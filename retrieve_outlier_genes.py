# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 08:50:24 2022

@author: Tyler
"""

from Bio import SeqIO
import pandas as pd

fasta_file = "mywa_2.1.fna" # Input fasta file
wanted_file = "gene_match_overlap.txt"
result_file = "overlap_gene_seqs.fa" # Output fasta file
outlier_type = "overlap"

genes = pd.read_table(wanted_file)

fasta_sequences = SeqIO.to_dict(SeqIO.parse(open(fasta_file),'fasta'))

#print(fasta_sequences.get("CM027508.1")[32953655:32978565].format("fasta"))

with open(result_file, "w") as f:
    for i in range(len(genes["SNP"])):
        #print(i)
        #print(genes.loc[i, "SNP"])
        #print(genes.loc[i, "chromosome"])
        #print(genes.loc[i, "gene"])
        #print(genes.loc[i, "start"] - 1)
        #print(genes.loc[i, "end"] - 1)
        g_snp = genes.loc[i, "SNP"]
        g_chrom = genes.loc[i, "chromosome"]
        g_gene = genes.loc[i, "gene"]
        g_start = genes.loc[i, "start"] - 1 
        g_end = genes.loc[i, "end"] - 1
        record = fasta_sequences.get(g_chrom)[g_start:g_end]
        record.id = outlier_type + ' ' + g_snp + ' ' + g_gene + ':'
        print(record.id)
        SeqIO.write([record], f, "fasta")  