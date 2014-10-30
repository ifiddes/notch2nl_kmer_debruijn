import networkx as nx
from Bio import SeqIO


class DeBruijnGraph(nx.DiGraph):
    def __init__(self, kmer_size):
        nx.DiGraph.__init__(self)
        self.kmer_size = kmer_size
    def add_sequences(self, seqRecord):
        k = self.kmer_size-1 #node is k-1 mers
        for i in xrange(len(seqRecord)-self.kmer_size):
            record = seqRecord[i:i+k]
            record.id = str(i)
            self.add_node(record)
    def build_edges(self):
        """
        Build all the edges of the graph. The nodes of the graph
       are k-1-mers. An edge between k-1-mers A and B is drawn if
       the end A is equal to the beginning of B
        """
        k = self.kmer_size-1
        for n1 in self.nodes_iter():
            for n2 in self.nodes_iter():
                if n1[1:k] == n2[0:k-1]:
                    self.add_edge(n1,n2)

    def prune_graph(self):
        """
        For each node, if has more than one outgoing edge remove all outgoing edges.
        Do the same for incoming edges.
        """
