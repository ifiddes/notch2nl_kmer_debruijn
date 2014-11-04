import networkx as nx
from Bio import SeqIO
import logging
#I wanted to have this class inherit the nx.DiGraph() directly
#but if I do, then finding subgraphs doesn't work correctly because I
#want to have kmer_size be initialized on instantiation, and then
#the copy function fails. SO has failed to help me so far.


class DeBruijnGraph(object):
    """
    Represents a DeBruijnGraph using networkx. Sequences are stored on each node as a
    k-1mer with edges pointing towards the next k-1mer in the original sequence.

    When initialized, a kmer size must be provided.

    To add sequences to this graph, pass Biopython SeqRecords to add_sequences().
    Once you have loaded sequences, prune the graph using prune_graph().

    prune_graph() removes all edges that fit into the rules below. This creates linear
    weakly connected components which represent continuous segments of sequence
    that represent one or more paralogs and can be used to infer copy number.

    1) every node with more than 1 outgoing edge has all outgoing edges removed
    2) every node with more than 1 incoming edge has all incoming edges removed
    3) every edge connecting two nodes with a different # of input sequences is removed
    This by definition removes all cycles as well as all strong connected components.

    """
    def __init__(self, kmer_size):
        logging.info("Initializing a DeBruijnGraph")

        self.kmer_size = kmer_size
        self.G = nx.DiGraph()
        self.has_sequences = False
        self.is_pruned = False


    def add_sequences(self, seqRecord):
        """
        Adds k1mers to the graph from a Biopython seqRecord. Edges are built as
        the k1mers are loaded.

        """
        k = self.kmer_size - 1
        name = seqRecord.name

        for i in xrange(len(seqRecord)-k):
            #left and right k-1mers
            km1L, km1R = str(seqRecord.seq[i:i+k]), str(seqRecord.seq[i+1:i+k+1])

            if self.G.has_node(km1L) is not True:
                self.G.add_node(km1L, pos=[i], source=[name], count=1)
            else:
                self.G.node[km1L]["pos"].append(i)
                self.G.node[km1L]["source"].append(name)
                self.G.node[km1L]["count"] += 1
            
            if self.G.has_node(km1R) is not True:
                self.G.add_node(km1R, pos=[], source=[], count=0)

            self.G.add_edge(km1L, km1R)

        #need to count the last kmer also
        self.G.node[km1R]["pos"].append(i+1)
        self.G.node[km1R]["source"].append(name)
        self.G.node[km1R]["count"] += 1

        self.has_sequences = True


    def prune_graph(self):
        """
        For each node, if has more than one outgoing edge remove all outgoing edges.
        Do the same for incoming edges. Also, remove all edges between nodes with
        different counts.

        """
        for n in self.G.nodes_iter():
            if len(self.G.in_edges(n)) > 1:
                for n1, n2 in self.G.in_edges(n):
                    self.G.remove_edge(n1, n2)

            if len(self.G.out_edges(n)) > 1:
                for n1, n2 in self.G.out_edges(n):
                    self.G.remove_edge(n1, n2)

            for n1, n2 in self.G.out_edges(n):
                if self.G.node[n1]["count"] != self.G.node[n2]["count"]:
                    self.G.remove_edge(n1, n2)

        self.is_pruned = True


    def weakly_connected_subgraphs(self):
        """
        Yields weakly connected subgraphs and their topolgical sort.

        """
        for subgraph in nx.weakly_connected_component_subgraphs(self.G):
            yield (subgraph, nx.topological_sort(subgraph))