import argparse, sys, os, logging
from itertools import izip
import cPickle as pickle

from src.deBruijnGraph import DeBruijnGraph
from Bio import SeqIO

"""
Builds a kmer DeBruijnGraph with nodes flagged for kmers present in genome.
Serializes this to disk for use by DeBruijnKmerILP.py

"""

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", "-r", type=str, required=True, 
        help="Reference fasta file")
    parser.add_argument("--out", "-o", type=argparse.FileType("wb"), 
        help="File to write pickled DeBruijnGraph to. Default is 'graphs/dbg.pickle'")
    parser.add_argument("--kmer_size", "-k", type=int, default=50, 
        help="kmer size. Default=50")
    parser.add_argument("--genome_counts", "-g", type=str, 
        help="Jellyfish kmer count fasta over genome sequences NOT containing region of interest. \
        Counts should be of k-1mers (49mer).")
    return parser.parse_args()


def parse_jellyfish_counts(file_path):
    with open(file_path) as f:
        for count, seq in izip(*[f]*2):
            yield seq.rstrip()


def main(args):
    args = parse_args(args)

    logging.basicConfig(filename="build_log.txt", format='%(asctime)-4s %(levelname)-6s %(message)s',
                    datefmt='%m-%d %H:%M', level=logging.DEBUG, filemode="w")
    
    logging.info("Building DeBruijnGraph.")
    G = DeBruijnGraph(args.kmer_size)
    paralogs = []
    
    for seqRecord in SeqIO.parse(args.reference, "fasta"):
        G.add_sequences(seqRecord)
        paralogs.append(seqRecord.name)
    
    G.prune_graph()

    logging.info("Reading genome kmer counts.")
    if args.genome_counts is not None:
        G.flag_nodes(parse_jellyfish_counts(args.genome_counts))

    logging.info("Dumping graph to disk.")
    pickle.dump(G, args.out)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
