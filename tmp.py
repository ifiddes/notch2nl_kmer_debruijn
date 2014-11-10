from src.kmerModel import KmerModel
from src.deBruijnGraph import DeBruijnGraph

from Bio import SeqIO
from itertools import izip
import cPickle as pickle
import argparse, sys, os, logging

G = pickle.load(open("graphs/dbg.pickle","rb"))

with open("test_data/snyder_test.fa") as f:
    data_counts = {seq.rstrip() : int(count.rstrip().lstrip(">")) for count, seq in izip(*[f]*2)}

outf = open("snyder_filtered.fa","w")
for kmer, count in data_counts.iteritems():
    if kmer in G.nodes():
        outf.write(">{}\n{}\n".format(count,kmer))

outf.close()
