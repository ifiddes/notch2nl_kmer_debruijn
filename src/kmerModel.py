from abstractIlpSolving import SequenceGraphLpProblem
from collections import Counter
from itertools import izip
import pulp, logging

class Block(object):
    """
    Represents one block (one weakly connected component)
    Each block may have anywhere from 1 to n paralogs represented in
    it. Initialized by passing in all of the nodes from one DeBruijnGraph
    WCC. Stores a mapping between each paralog this WCC represents
    to the start and stop positions in the source sequence.
    Also stores all of the kmers in this block as a set.

    """
    def __init__(self, subgraph, topo_sorted, default_ploidy=2):
        #a set of all kmers represented by this block
        self.kmers = set()
        for kmer in topo_sorted:
            #blocks do not contain kmers that are elsewhere in the genome
            if 'bad' in subgraph.node[kmer]:
                self.kmers.add(kmer)

        #number of bases in this window
        if len(self.kmers) > 0:
            self.size = len(self.kmers) + len(kmer) - 1
        else:
            self.size = 0

        #a mapping of paralog, position pairs to variables
        #one variable for each instance of a input sequence
        self.variables = {}

        #since each node has the same set of sequences, and
        #we have a topological sort, we can pull down the positions
        #of the first and last node
        
        start_kmer, stop_kmer = topo_sorted[0], topo_sorted[-1]
        start_node, stop_node = subgraph.node[start_kmer], subgraph.node[stop_kmer]
        
        #build variables for each instance
        for para, start, stop in izip(start_node["source"], start_node["pos"], stop_node["pos"]):
            stop = stop + len(kmer) #off by one?
            #restrict the variable to be an integer with default_ploidy
            if self.get_size() != 0:
                self.variables[(para, start, stop)] = pulp.LpVariable("{}:{}-{}".format(
                    para, start, stop), default_ploidy, cat="Integer")
            else:
                self.variables[(para, start, stop)] = None

    def kmer_iter(self):
        """iterator over all kmers in this block"""
        for kmer in self.kmers:
            yield kmer

    def variable_iter(self):
        """iterator over variables and related values in this block"""
        for (para, start, stop), variable in self.variables.iteritems():
            yield para, start, stop, variable

    def get_variables(self):
        """returns all LP variables in this block"""
        return self.variables.values()

    def get_size(self):
        """returns size of block"""
        return self.size

    def get_kmers(self):
        """returns set of all kmers in block"""
        return self.kmers


class KmerModel(SequenceGraphLpProblem):
    """
    Represents a kmer DeBruijnGraph model to infer copy number of highly
    similar paralogs using ILP.

    Each weakly connected component of the graph represents one 'block'
    which will be assigned one LP variable for each source sequence
    within the WCC (each node all comes from the same sequence(s) by the
    pruning method used).
    
    All constraints and penalty terms are automatically added to the
    SequenceGraphLpProblem to which the Model is attached (specified on
    construction).
    
    """
    def __init__(self, paralogs):
        logging.debug("Initializing KmerModel.")
        SequenceGraphLpProblem.__init__(self)
        self.blocks = []
        self.block_map = { x : [] for x in paralogs }
        self.built = False
        self.has_data = False

    def has_data(self):
        return self.has_data

    def is_built(self):
        return self.built


    def build_blocks(self, DeBruijnGraph, breakpoint_penalty=25):
        """
        Builds a ILP kmer model. Input:

        DeBruijnGraph, which is a networkx DeBruijnGraph built over the genome region of interest.

        breakpoint_penalty determines the penalty that ties instances together

        """
        #make sure the graph has been initialized and pruned
        assert DeBruijnGraph.is_pruned and DeBruijnGraph.has_sequences
        logging.debug("Building blocks in KmerModel.")

        #build the blocks, don't tie them together yet
        for subgraph, topo_sorted in DeBruijnGraph.weakly_connected_subgraphs():
            b = Block(subgraph, topo_sorted)
            self.blocks.append(b)

        logging.debug("Blocks built.")

        for block in self.blocks:
            for para, start, stop, variable in block.variable_iter():
                self.block_map[para].append([start, stop, variable])

        logging.debug("block_map built.")

        #now sort these maps and start tying variables together
        for para in self.block_map:
            self.block_map[para] = sorted(self.block_map[para], key = lambda x: x[0])

            #check for off by 1 here?
            #filter out all blocks without variables (no kmers)
            variables = [v for s, st, v in self.block_map[para] if v is not None]

            for i in xrange(1, len(variables)):
                var_a, var_b = variables[i-1], variables[i]
                self.constrain_approximately_equal(var_a, var_b, breakpoint_penalty)

        logging.debug("Block variables tied together; block_map sorted.")
        self.built = True


    def introduce_data(self, kmerCounts, data_penalty=50):
        """
        Introduces data to this ILP kmer model. For this, the input is assumed to be a dict 
        representing the results of kmer counting a WGS dataset (format seq:count)

        data_penalty represents the penalty in the ILP model for the copy number of each instance
        to deviate from the data.

        coverage represents the best guess of per-base coverage across the genome for this dataset.
        This will be used to normalize the counts to account for variable block sizes. Getting this
        number correct is kind of important.

        """
        logging.debug("Starting to introduce {} kmers to model.".format(len(kmerCounts)))
        i = 0

        for block in self.blocks:
            i += 1
            if i % 100 == 0:
                logging.debug("Examined {} windows".format(i))

            if block.get_size() > 0:
                count = sum(kmerCounts.get(x, 0) for x in block.get_kmers())
                
                #var_b is the sum of all LP variables in the window
                window_lp_sum = self.penalties.get_variable()
                
                #constrain var_b as equal to the sum of LP variables (only way to constrain sum)
                self.problem += window_lp_sum == sum(block.get_variables())
                
                #constrain approximately equal subject to a data_penalty
                self.constrain_approximately_equal(count, window_lp_sum, data_penalty * block.get_size())

        self.has_data = True

    def report_copy_number(self):
        """
        Reports copy number from solved ILP problem. Loops over the block_map class member
        and reports the linear copy number calls for each paralog.
        """
        assert self.is_solved
        logging.debug("Reporting copy number from solved problem")

        copy_map = { x : [] for x in self.block_map.keys()}

        #self.block_map is sorted by start position, so we can just loop through it
        for para in self.block_map:
            for start, stop, var in self.block_map[para]:
                if var is not None:
                    c = pulp.value(var)
                    copy_map[para].append([start, c])

        return copy_map