# mutation-detection
- To be able to use this script, you should first align the reads required to be checked for mutations against wildtype polypeptides using BLASTX:

    blastx -query test.fasta -db wt_polypeptide -num_descriptions 1 -num_alignments 1 -out test.vs.wt_polypeptide.blastx.txt

- You should also prepare a mutation list (see example_list.txt).
- The mutation list is a tab-separated text file with first column containing protein idntifier, second column containing amino acid position of the mutation and the third column contains the one-letter abbreviation of the mutated amino acid (NOT the wildtype).
- Use the mutation list and the plain BLAST report as inputs for the script.
- Usage: ./mutation-detection-script.pl [input mutation list file] [input blast text report]
- The script reports reads containing mutations included in the list to standard output (see example_output.txt).
