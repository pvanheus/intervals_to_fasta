# intervals_to_fasta
Given a GFF3 description of genes and CDS on a sequence and a FASTA file of the sequence, output the intergenic intervals as FASTA format file. Similar to @Friedberg-Lab's intergenic.py

Usage:

usage: intervals_to_fasta.py [-h] [--pad_start PAD_START] [--pad_end PAD_END]
                             fasta_input_file gff_input_file [output_file]

**NOTE:** requires pvh's gff_utils module