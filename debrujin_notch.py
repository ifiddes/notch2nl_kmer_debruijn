from collections import defaultdict
import string, sys, os, argparse
from Bio import SeqIO
"""
Implementation of a Naive DeBrujin graph with modifications
to retain information about which input sequence provided the kmer
as well as which position it is in that sequence.
Code adapted from https://github.com/rchikhi/python-bcalm
"""

revcomp_trans=string.maketrans('ATGC', 'TACG')

class Node(object):
    """
    Represents one Node in the DeBrujin graph.
    Each node stores the kmer sequence, the position in the
    reference sequence it came from, as well as the name
    of the sequence it came from as it is in the fasta.
    """
    __slots__ = ['seq','start','source'] #save memory
    def __init__(self, seq, pos, source):
        self.seq = seq
        self.start = pos
        self.source = source
    def __len__(self):
        """overload so len() works"""
        return len(self.seq)
    def __getitem__(self, key):
        """overload so slicing works"""
        if isinstance(key, int) or isinstance(key, slice):
            return self.seq[key]
        else:
            raise AssertionError("tried to slice a Node with something that is not a int or a slice")
    def rc(self):
        """reverse complement the kmer"""
        return string.translate(self.seq, revcomp_trans)[::-1]

class Graph:
    """
    Represents a DeBrujin graph (although it is not a DeBrujin graph until
    the debrujin() function has been run).
    Initialize then use loadkmers() to load kmers from a fasta, then run
    debrujin() to make the edges.
    """
    def __init__(self, k):
        self.k = k
        self.nodes = dict()
        self.last_node_index = 0
        self.map = defaultdict(list) # last k-1-mers
        self.maprev = defaultdict(list) # last k-1-mers, reverse complemented
        self.neighbor = defaultdict(list)
        
    def addvertex(self, vertex):
        """Add a vertex to the graph"""
        key = vertex[:self.k-1]
        keyrc = vertex.rc()[:self.k-1]
        idx = self.last_node_index 
        self.nodes[idx] = vertex
        self.map[key] += [idx]
        self.maprev[keyrc] += [idx]
        self.last_node_index += 1

    def loadkmers(self, name, iterable):
        """Used by the main and kmerize() to load kmers as nodes"""
        nb_nodes, nb_nt = 0, 0
        for seq, pos in iterable: #loop over kmers
            n = Node(seq, pos, name)
            self.addvertex(n)

    def add_edge(self, i1, i2, label):
        """Add a edge between two vertices"""
        self.neighbor[i1] += [(i2, label[0])]
        self.neighbor[i2] += [(i1, label[1])]

    def debruijn(self):
        """Finds edges between nodes in the graph"""
        for i, node in self.nodes.items():
            key = node[-(self.k-1):]
            for other_idx in self.map[key]:
                self.add_edge(i, other_idx, 'OI') # I'm using the In/Out notation, this is a >--> edge
            for other_idx in self.maprev[key]:
                # prevent adding same edge twice
                if i < other_idx:
                    self.add_edge(i, other_idx, 'OO') # >---< edge
            keyrc = node.rc()[-(self.k-1):]
            for other_idx in self.map[keyrc]:
                if i < other_idx:
                    self.add_edge(i, other_idx, 'II') # <---> edge
        del self.map
        del self.maprev



def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--files", "-f", nargs="+", type=str, required=True, help="1 or more fasta files")
    parser.add_argument("--kmer_size", "-k", type=int, default=50, help="kmer size. Default=50")
    return parser.parse_args()


def kmerize(record, k):
    """Generator to yield (kmer, position) tuples over a biopython seqRecord"""
    s = str(record.seq).upper()
    for i in xrange(len(record)-k):
        yield s[i:i+k], i

def main(args):
    args = parse_args(args)
    G = Graph(args.kmer_size)
    for f in args.files:
        for record in SeqIO.parse(f,"fasta"):
            G.loadkmers(record.id, kmerize(record, args.kmer_size))
    G.debruijn()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
