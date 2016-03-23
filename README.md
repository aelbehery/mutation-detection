# mutation-detection
- To be able to use this script, you should first align the reads required to be checked for mutations against wildtype polypeptides using BLASTX:

    blastx -query test.fasta -db wt_polypeptide -num_descriptions 1 -num_alignments 1 -out test.vs.wt_polypeptide.blastx.txt

- You should also prepare a mutation list. This list is a tab-separated text file with the first column containing protein idntifier, the second one containing amino acid position of the mutation and the third one containing the one-letter abbreviation of the mutated amino acid NOT the wildtype (see example_list.txt).
- Use the mutation list and the plain BLAST report as inputs for the script.
- Usage: ./mutation-detection-script.pl [input mutation list file] [input blast text report]
- The script reports reads containing mutations included in the list to standard output (see example_output.txt).
