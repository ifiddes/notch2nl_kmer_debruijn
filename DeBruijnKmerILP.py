#!/usr/bin/env python2.7

#from jobTree.scriptTree.target import Target 
#from jobTree.scriptTree.stack import Stack 

from src.kmerModel import KmerModel
from src.deBruijnGraph import DeBruijnGraph

from Bio import SeqIO
from itertools import izip
import cPickle as pickle
import argparse, sys, os, logging, string, gzip


def parse_args(args):
    parser = argparse.ArgumentParser()
    cwd = os.getcwd()
    parser.add_argument("--graph", "-g", type=argparse.FileType("rb"), 
        help="Pickled DeBruijnGraph to load. Default is 'graphs/dbg.pickle'",
        default=os.path.join(cwd,"notch2nl_kmer_debruijn/graphs/dbg.pickle"))
    parser.add_argument("--sample_counts", "-c", type=str, required=True, 
        help="Jellyfish k-1mer (49bp) counts for this sample.")
    parser.add_argument("--normalizing_counts", "-n", type=int, required=True,
        help="Samtools view -c output over normalizing region for this sample.")
    parser.add_argument("--normalizing_size", "-s", type=int, required=True, default=1408091,
        help="Normalizing size that --normalizing_counts was drawn from.")
    parser.add_argument("--breakpoint_penalty", "-b", type=float, 
        help="Breakpoint penalty for ILP.", default=2500)
    parser.add_argument("--data_penalty", "-d", type=float, 
        help="Data penalty for ILP.", default=500)
    return parser.parse_args()


def write_result(copy_map):
    """tmp debugging write"""
    for para in copy_map:
        outf = open(os.path.join("output", "tmp_{}.txt".format(para)), "w")
        for x in copy_map[para]:
            outf.write("\t".join(map(str,x)) + "\n")
        outf.close()

def main(args):
    args = parse_args(args)

    logging.basicConfig(filename="log.txt", format='%(asctime)-4s %(levelname)-6s %(message)s',
                    datefmt='%m-%d %H:%M', level=logging.DEBUG, filemode="w")
    
    logging.info("Loading DeBruijnGraph.")
    G = pickle.load(args.graph)

    logging.info("Building data counts.")
    with gzip.open(args.sample_counts, "rb") as f:
        #using seq.translate has been shown to be faster than using any other means
        #of removing characters from a string
        rm = ">\n"
        #kmers = G.kmers()
        data_counts = {seq.translate(None, rm) : int(count.translate(None, rm)) for
            count, seq in izip(*[f]*2) if seq.translate(None, rm) in G.kmers}

    coverage = 1.0 * normalizing_counts / normalizing_size

    #if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
    #    pickle.dump(data_counts,open("tmp_data_counts.pickle","wb"))

    logging.info("Initializing kmer model.")
    P = KmerModel(G.paralogs, coverage)
    P.build_blocks(G, args.breakpoint_penalty)

    logging.info("Introducing kmer count data to model.")
    P.introduce_data(data_counts, args.data_penalty)

    logging.info("Solving ILP model.")
    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        P.solve(save="test.lp")
    else:
        P.solve()

    logging.info("Solved! Getting copy number.")
    copy_map = P.report_copy_number()

    write_result(copy_map)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
