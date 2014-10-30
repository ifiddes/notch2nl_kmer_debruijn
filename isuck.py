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

def kmerize(record, k):
    """Generator to yield (kmer, position) tuples over a biopython seqRecord"""
    s = str(record.seq).upper()
    for i in xrange(len(record)-k):
        yield s[i:i+k], i


from debruijn_notch import Graph

G = Graph(50)


for record in SeqIO.parse("notch2.fasta","fasta"):
   G.load_kmers(record.id, kmerize(record,50))

G.debruijn()
G.prune_graph()

