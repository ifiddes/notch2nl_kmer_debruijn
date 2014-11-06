#!/usr/bin/env python2.7

#from jobTree.scriptTree.target import Target 
#from jobTree.scriptTree.stack import Stack 

from src.kmerModel import KmerModel
from src.deBruijnGraph import DeBruijnGraph

from Bio import SeqIO
from itertools import izip
import argparse, sys, os, logging


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", "-r", type=str, required=True, 
        help="Reference fasta file")
    parser.add_argument("--kmer_filter", "-kf", type=str, 
        help="Optional jellyfish fasta file of kmers NOT in reference region.")
    parser.add_argument("--kmer_size", "-k", type=int, default=50, 
        help="kmer size. Default=50")
    parser.add_argument("--sample_counts", "-c", type=str, required=True, 
        help="Jellyfish kmer counts for sample whose CNV we are interested in.")
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
        outf = open(os.path.join(output, "tmp_{}.txt".format(para)), "w")
        for x in copy_map[para]:
            outf.write("\t".join(map(str,x)) + "\n")
        outf.close()

def main(args):
    args = parse_args(args)

    logging.basicConfig(filename="log.txt", format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M', level=logging.DEBUG, filemode="w")
    logging.info("Building DeBruijnGraph (main)")

    #build the DeBruijnGraph
    G = DeBruijnGraph(args.kmer_size)
    paralogs = []
    for seqRecord in SeqIO.parse(args.reference, "fasta"):
        G.add_sequences(seqRecord)
        paralogs.append(seqRecord.name)
    G.prune_graph()

    #build the kmer_filter set if provided (likely to be VERY slow with HUGE RAM requirements)
    if args.kmer_filter is not None:
        logging.info("Building kmer_filter (main)")
        kmer_filter = set()
        #not using biopython for speed
        with open(args.kmer_filter) as f:
            add = kmer_filter.add
            for length, seq in izip(*[f]*2):
                add(seq.rstrip())
    else:
        kmer_filter = None

    logging.info("Building data counts (main)")
    #build a dict of kmer counts seen the WGS extracted region + unmapped
    #also not using biopython for speed
    with open(args.sample_counts) as f:
        data_counts = {seq.rstrip() : int(count.rstrip().lstrip(">")) for count, seq in izip(*[f]*2)}

    logging.info("Initializing kmer model (main)")
    #now we initialize and populate the kmer model
    P = KmerModel(paralogs)
    P.build_blocks(G, args.breakpoint_penalty)

    #introduce the data and solve
    P.introduce_data(data_counts, kmer_filter, args.coverage, args.data_penalty)
    logging.info("Solving... (main)")
    P.solve(save="tmp_ilp.txt")

    logging.info("Solved! Getting copy number (main)")
    copy_map = P.report_copy_number()

    write_result(copy_map)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
