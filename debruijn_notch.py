from collections import defaultdict, Counter
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
        self.nodes = {}
        self.last_node_index = 0
        self.map = defaultdict(list) # last k-1-mers
        self.maprev = defaultdict(list) # last k-1-mers, reverse complemented
        self.edges = defaultdict(list)
        
    def add_vertex(self, vertex):
        """Add a vertex to the graph"""
        key = vertex[:self.k-1]
        keyrc = vertex.rc()[:self.k-1]
        idx = self.last_node_index 
        self.nodes[idx] = vertex
        self.map[key] += [idx]
        self.maprev[keyrc] += [idx]
        self.last_node_index += 1

    def load_kmers(self, name, iterable):
        """Used by the main and kmerize() to load kmers as nodes"""
        for seq, pos in iterable: #loop over kmers
            n = Node(seq, pos, name)
            self.add_vertex(n)

    def add_edge(self, i1, i2, label):
        """Add a edge between two vertices"""
        self.edges[i1] += [(i2, label[0])]
        self.edges[i2] += [(i1, label[1])]

    def debruijn(self):
        """Finds edges between nodes in the graph"""
        for i, node in self.nodes.iteritems():
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

    def prune_graph(self):
        """
        For each node, if has more than one outgoing edge remove all outgoing edges.
        Do the same for incoming edges.
        """
        for i, node in self.nodes.iteritems():
            directions = Counter([d for p, d in self.edges[i]])
            if directions["I"] > 1 and directions["O"] > 1:
                self.edges[i] = []
            elif directions["I"] > 1:
                self.edges[i] = [x for x in self.edges[i] if x[1] != "I"]
            elif directions["O"] > 1:
                self.edges[i] = [x for x in self.edges[i] if x[1] != "O"]

    def dfs(self):
        """
        Find all (weakly) connected componnents in the graph using depth-first search
        """
        


    def strongly_connected_components(self):
        """
        Find the strongly connected components in a graph using
        Tarjan's algorithm.
        Adapted from http://www.logarithmic.net/pfh-files/blog/01208083168/sort.py
        """
        result = []
        stack = []
        low = {}

        def visit(node):
            if node in low: return

            num = len(low)
            low[node] = num
            stack_pos = len(stack)
            stack.append(node)

            for successor in self.edges[node]:
                if len(self.edges[successor][1]) > 0:
                    if self.edges[successor][1] == "O":
                        visit(successor)
                        low[node] = min(low[node], low[successor])

            if num == low[node]:
                component = tuple(stack[stack_pos:])
                del stack[stack_pos:]
                result.append(component)
                for item in component:
                    low[item] = len(self.nodes)

        for node in self.nodes:
            visit(node)
        return result


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
            G.load_kmers(record.id, kmerize(record, args.kmer_size))
    G.debruijn()
    G.prune_graph()
    SCCs = G.strongly_connected_components()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
