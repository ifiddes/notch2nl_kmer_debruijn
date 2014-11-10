#!/usr/bin/env python2.7

#from jobTree.scriptTree.target import Target 
#from jobTree.scriptTree.stack import Stack 

from src.kmerModel import KmerModel
from src.deBruijnGraph import DeBruijnGraph

from Bio import SeqIO
from itertools import izip
import cPickle as pickle
import argparse, sys, os, logging


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", "-g", type=argparse.FileType("rb"), 
        help="Pickled DeBruijnGraph to load. Default is 'graphs/dbg.pickle'",
        default="graphs/dbg.pickle")
    parser.add_argument("--sample_counts", "-c", type=str, required=True, 
        help="Jellyfish k-1mer (49bp) counts for sample whose CNV we are interested in.")
    parser.add_argument("--breakpoint_penalty", "-b", type=float, 
        help="Breakpoint penalty for ILP.", default=50)
    parser.add_argument("--data_penalty", "-d", type=float, 
        help="Data penalty for ILP.", default=5)
    parser.add_argument("--coverage", "-cov", type=float, 
        help="Expected per-base coverage for this WGS.", default=30)
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
    with open(args.sample_counts) as f:
        data_counts = {seq.rstrip() : int(count.rstrip().lstrip(">")) for count, seq in izip(*[f]*2)}

    logging.info("Initializing kmer model.")
    P = KmerModel(paralogs)
    P.build_blocks(G, args.breakpoint_penalty)

    logging.info("Introducing kmer count data to model.")
    P.introduce_data(data_counts, args.coverage, args.data_penalty)

    logging.info("Solving ILP model.")
    P.solve(save="test.lp")

    logging.info("Solved! Getting copy number.")
    copy_map = P.report_copy_number()

    write_result(copy_map)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
