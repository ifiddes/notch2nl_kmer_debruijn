notch2nl_kmer_debruijn
======================
Code for doing debruijn graphs for notch2nl. Will have ILP code too eventually.

Right now, only `debruijn_network_not_oo.py` works. As it stands, this script will take a fasta file and a optional name and kmer size argument (default kmer size is 50, default file name header is component) and will build the graph, prune the edges, and print out a fasta file containing the sequence of each weakly connected component. It will also write a text file with the sizes of each of these and make a histogram. Finally, if you are on the hive filesystem and have `bwa` in your path, it will align these back to the CHM1 assembly and provide a BED file of these alignments. If you do not have BWA it will crash.

If you want to run a test version of this that will run on both all four notch2nl paralogs as well as just notch2, use run.sh.
