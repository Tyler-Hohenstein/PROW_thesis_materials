# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 17:32:09 2022

@author: Tyler
"""

from Bio import SeqIO

fasta_file = "reference.fasta" # Input fasta file
wanted_file = "overlap.outlier.contigs.txt" # Input interesting sequence IDs, one per line
result_file = "overlap.outlier.contigs.fa" # Output fasta file

wanted = set()
with open(wanted_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], f, "fasta")