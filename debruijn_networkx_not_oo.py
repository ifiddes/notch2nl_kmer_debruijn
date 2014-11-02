import networkx as nx
from Bio import SeqIO
from collections import defaultdict, Counter
import argparse, sys, os

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--files", "-f", nargs="+", type=str, required=True, help="1 or more fasta files")
    parser.add_argument("--kmer_size", "-k", type=int, default=50, help="kmer size. Default=50")
    parser.add_argument("--name", type=str, default="components", help="debugging file name header")
    return parser.parse_args()


def add_sequences(G, seqRecord, kmer_size):
    k = kmer_size - 1
    name = seqRecord.name
    for i in xrange(len(seqRecord)-k):
        #left and right k-1mers
        km1L, km1R = str(seqRecord.seq[i:i+k]), str(seqRecord.seq[i+1:i+k+1])
        if G.has_node(km1L) is not True:
            G.add_node(km1L, pos=[i], source=[name], count=1)
        else:
            G.node[km1L]["pos"].append(i); G.node[km1L]["source"].append(name); G.node[km1L]["count"] += 1
        if G.has_node(km1R) is not True:
            G.add_node(km1R, pos=[], source=[], count=0)
            #we initialize the right k1mer as a node to have a edge
            #but we don't give it a count/source/pos because then it would be double counted
            #on the next iteration
        G.add_edge(km1L, km1R)
    #need to count the last kmer also
    G.node[km1R]["pos"].append(i+1); G.node[km1R]["source"].append(name); G.node[km1R]["count"] += 1


def prune_graph(G):
    """
    For each node, if has more than one outgoing edge remove all outgoing edges.
    Do the same for incoming edges. Also, remove all edges between nodes with
    different counts.
    """
    for n in G.nodes_iter():
        if len(G.in_edges(n)) > 1:
            for n1, n2 in G.in_edges(n):
                G.remove_edge(n1, n2)
        if len(G.out_edges(n)) > 1:
            for n1, n2 in G.out_edges(n):
                G.remove_edge(n1, n2)
        for n1, n2 in G.out_edges(n):
            if G.node[n1]["count"] != G.node[n2]["count"]:
                G.remove_edge(n1, n2)


def write_weak_subgraphs_to_fasta_and_sizes(G, name):
    """temp debugging code"""
    outf = open("{}_components.fasta".format(name),"w")
    outsizes = open("{}_component_sizes.txt".format(name),"w")
    i=0
    tmp = []

    for subgraph in nx.weakly_connected_component_subgraphs(G):
        ordered = nx.topological_sort(subgraph)
        seq = [ordered[0]]
        names = ",".join(["_".join(map(str,x)) for x in Counter(G.node[ordered[0]]["source"]).items()])
        for x in ordered[1:]:
            seq.append(x[-1])
        tmp.append("".join(seq))
        outf.write(">{}\n{}\n".format(str(i)+"|"+names,"".join(seq)))
        outsizes.write("{}\t{}\n".format(str(i)+"|"+names,len(seq)))
        i+=1

    outf.close()
    outsizes.close()

def align_to_chm1(name):
    os.system("bwa mem /cluster/home/ifiddes/chr1_fastas/chm1_chr1.fa components_notch2_all_paralogs.fasta | samtools view -bS - | bedtools bamtobed > {}.bed".format(name))


def main(args):
    args = parse_args(args)
    G = nx.DiGraph()

    for f in args.files:
        for record in SeqIO.parse(f,"fasta"):
            add_sequences(G, record, args.kmer_size)

    prune_graph(G)

    write_weak_subgraphs_to_fasta_and_sizes(G, args.name)
    align_to_chm1(args.name)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
