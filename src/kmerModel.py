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
        #save the count for debugging, after we load data
        self.count = 0
        for kmer in topo_sorted:
            #blocks do not contain kmers that are elsewhere in the genome
            #if 'bad' not in subgraph.node[kmer]:
            self.kmers.add(kmer)

        #size of this window
        self.size = len(self.kmers) + len(kmer)

        #a mapping of paralog, position pairs to variables
        #one variable for each instance of a input sequence
        self.variables = {}

        #since each node has the same set of sequences, and
        #we have a topological sort, we can pull down the positions
        #of the first node only
        start_node = subgraph.node[topo_sorted[0]]
        
        #build variables for each instance
        for para, start in izip(start_node["source"], start_node["pos"]):
            #restrict the variable to be an integer with default_ploidy
            if len(self.kmers) > 0:
                self.variables[(para, start)] = pulp.LpVariable("{}:{}".format(
                    para, start), default_ploidy, cat="Integer")
            else:
                self.variables[(para, start)] = None

    def kmer_iter(self):
        """iterator over all kmers in this block"""
        for kmer in self.kmers:
            yield kmer

    def variable_iter(self):
        """iterator over variables and related values in this block"""
        for (para, start), variable in self.variables.iteritems():
            yield para, start, variable

    def get_variables(self):
        """returns all LP variables in this block"""
        return self.variables.values()

    def get_size(self):
        """returns size of block"""
        return self.size

    def get_kmers(self):
        """returns set of all kmers in block"""
        return self.kmers

    def get_count(self):
        """returns counts of data seen in this block"""
        return self.count


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


    def build_blocks(self, DeBruijnGraph, breakpoint_penalty=500):
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
            for para, start, variable in block.variable_iter():
                self.block_map[para].append([start, variable, block])

        logging.debug("block_map built.")

        #now sort these maps and start tying variables together
        for para in self.block_map:
            self.block_map[para] = sorted(self.block_map[para], key = lambda x: x[0])

            #filter out all blocks without variables (no kmers)
            #store block sizes for normalizing purposes
            variables = [[v, b.get_size()] for s, v, b in self.block_map[para] if v is not None]

            for i in xrange(1, len(variables)):
                var_a, var_b = variables[i-1][0], variables[i][0]
                a_norm, b_norm = variables[i-1][1], variables[i][1]
                self.constrain_approximately_equal(var_a / a_norm, 
                    var_b / b_norm, breakpoint_penalty)

        logging.debug("Block variables tied together; block_map sorted.")
        self.built = True


    def introduce_data(self, kmerCounts, data_penalty=50):
        """
        Introduces data to this ILP kmer model. For this, the input is assumed to be a dict 
        representing the results of kmer counting a WGS dataset (format seq:count)

        data_penalty represents the penalty in the ILP model for the copy number of each instance
        to deviate from the data, and is scaled to the size of the block.

        """
        logging.debug("Starting to introduce {} kmers to model.".format(len(kmerCounts)))

        for block in self.blocks:
            if len(block.get_kmers()) > 0:
                count = sum(kmerCounts.get(x, 0) for x in block.get_kmers())
                #constrain approximately equal subject to a data_penalty
                block_penalty = data_penalty * block.get_size()
                block.count = count
                self.constrain_approximately_equal(count, sum(block.get_variables(), block_penalty))

        self.has_data = True

    def report_copy_number(self):
        """
        Reports copy number from solved ILP problem. Loops over the block_map class member
        and reports the linear copy number calls for each paralog.
        """
        assert self.is_solved
        logging.debug("Reporting copy number from solved problem")

        copy_map = { x : [] for x in self.block_map.keys()}

        for para in self.block_map:
            for start, var, block in self.block_map[para]:
                size = block.get_size()
                count = block.get_count()
                num_vars = len(block.get_variables())
                if var is not None:
                    c = pulp.value(var)
                    copy_map[para].append([start, size-49, count, num_vars, var._LpElement__name]) 
                else:
                    copy_map[para].append([start, size-49, count, num_vars, None]) 

        return copy_map