import networkx as nx
from Bio import SeqIO
from collections import defaultdict
import argparse, sys, os
"""
This version doesn't work because I have not found an effective way to overload
nx.DiGraph.subgraph() such that nx.weakly_connected_component_subgraphs will work.
See: http://stackoverflow.com/questions/26681892/overloading-subgraph-so-that-derived-graph-objects-can-have-subgraphs-be-made-n/26690286#26690286

"""

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--files", "-f", nargs="+", type=str, required=True, help="1 or more fasta files")
    parser.add_argument("--kmer_size", "-k", type=int, default=50, help="kmer size. Default=50")
    parser.add_argument("--name", type=str, default="components", help="debugging file name header")
    return parser.parse_args()

class DeBruijnGraph(nx.DiGraph):
    """
    Represents a deBruijnGraph for inferring copy number of highly related paralogs.
    Initialize by passing seqRecords to add_sequences(), then prune the graph with prune_graph()
    to remove all edges between kmers with variable copy number as well as kmers with
    multiple incoming or outgoing edges.
    """
    def __init__(self, kmer_size):
        nx.DiGraph.__init__(self)
        self.kmer_size = kmer_size
    

    def add_sequences(self, seqRecord):
        k = self.kmer_size - 1
        name = seqRecord.name
        
        for i in xrange(len(seqRecord)-k):
            #left and right k-1mers
            km1L, km1R = str(seqRecord.seq[i:i+k]), str(seqRecord.seq[i+1:i+k+1])
            
            if self.has_node(km1L) is not True:
                self.add_node(km1L, pos=[i], source=[name], count=1)
            else:
                self.node[km1L]["pos"].append(i); self.node[km1L]["source"].append(name); self.node[km1L]["count"] += 1
            
            if self.has_node(km1R) is not True:
                self.add_node(km1R, pos=[], source=[], count=0)
            self.add_edge(km1L, km1R)
        
        #need to count the last kmer also
        self.node[km1R]["pos"].append(i+1); self.node[km1R]["source"].append(name); self.node[km1R]["count"] += 1
    

    def prune_graph(self):
        """
        For each node, if has more than one outgoing edge remove all outgoing edges.
        Do the same for incoming edges. Also, remove all edges between nodes with
        different counts.
        """
        for n in self.nodes_iter():
            if len(self.in_edges(n)) > 1:
                for n1, n2 in self.in_edges(n):
                    self.remove_edge(n1, n2)

            if len(self.out_edges(n)) > 1:
                for n1, n2 in self.out_edges(n):
                    self.remove_edge(n1, n2)

            for n1, n2 in self.out_edges(n):
                if self.node[n1]["count"] != self.node[n2]["count"]:
                    self.remove_edge(n1, n2)

    def subgraph(self, nbunch, copy=True):
        H = self.subgraph(self, nbunch, copy)
        H.kmer_size = self.kmer_size
        if copy:
            return DeBruijnGraph(H)
        else:
            return H


def write_graph_to_fasta_and_sizes(G, name):
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
    G = DeBruijnGraph(args.kmer_size)

    for f in args.files:
        for record in SeqIO.parse(f,"fasta"):
            G.add_sequences(record)
    G.prune_graph()
    weak_subgraphs = nx.weakly_connected_component_subgraphs(G)

    write_graph_to_fasta_and_sizes(G, args.name)
    align_to_chm1(args.name)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
