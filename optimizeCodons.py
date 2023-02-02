#!/usr/bin/env python
# coding: utf-8

# # optimizeCodons.py
# 
# Back-translate protein sequences by semi-randomized codon optimization strategy.
# 
# Inputs:
# - File, Protein multiple sequence alignment (MSA) (FASTA format), containing the target and template protein sequence(s)
# - File, Codon usage table (EMBOSS format) for organism with template protein
# - File, Nucleotide sequence (FASTA format) of template protein (including stop codon(s))
# - Python list, specifying alignment column intervals per protein in MSA
# - (if needed) regex pattern for target sequences
#         
# Method:
# - Compare target and template sequence
# - For sites where target and template are identical, assign template codon.
# - For sites where target and template differ, assign codon by random choice, weighted by codon usage table.
# - Compatible with alignments containing concatenated proteins, given user-specified alignment column intervals.
# - Adds start codon if target protein(s) missing N-terminal Met residue
# 
# Output:
# - Fasta file containing codon-optimized nucleotide sequences of target proteins
# 
# Last updated: 2023-01-26, Amanda K. Garcia (akgarcia3@wisc.edu)

# ### Setup

import re
import pandas as pd
import random
import os


# ### Import data

# Assign starting variables

msa_f = "align.fasta"
template_dna_f = "template_DNA.fasta"
codon_table_f = "codonusagetable.cusp"
template = "template_name"
# target_pattern = "" # Specify target sequence(s)
gene_delim = [0, 449, 1392, 2502] # Specify alignment column intervals defining proteins in concatenated alignment.
                                  # Format: [first column, last column of gene 1, last column of gene 2, ..., last column]
                                  # If single protein alignment, specify only first and last columns


# Import target and template protein sequences from MSA

with open(msa_f, "rt") as myfile: # Parse MSA
    msa_lines = myfile.readlines()

template_aa = "" # Template protein sequence
for line in msa_lines:
    if template in line:
        template_aa = re.sub("\n|\r", "", msa_lines[msa_lines.index(line) + 1])
        template_aa_split = [template_aa[i:j] for i,j in zip(gene_delim, gene_delim[1:]+[None])] # Add "Stop" characters, "*"
        template_aa = "*".join(template_aa_split)

target_msa = [] # List containing target protein name/sequence pair(s)
for line in msa_lines:
    if ">" in line: # Replace ">" with different target pattern here if needed
        seq_name = re.sub(">|\n|\r", "", line)
        seq = re.sub("\n|\r", "", msa_lines[msa_lines.index(line) + 1])
        seq_split = [seq[i:j] for i,j in zip(gene_delim, gene_delim[1:]+[None])] # Add "Stop" characters, "*"
        seq = "*".join(seq_split)
        target_msa.append(tuple([seq_name, seq]))


# Import template DNA sequence

template_dna = "" # Template DNA sequence
with open(template_dna_f, "rt") as myfile:
    template_dna = myfile.readlines()[1]
    template_dna = re.sub("\n|\r","",template_dna)


# Import codon usage table

df = pd.read_csv(codon_table_f, delim_whitespace=True)

aa_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]

codon_dict = {} # dictionary containing codon/frequency list pairs for each amino acid

codon_lists = []

for res in aa_list:
    codon_list = df.loc[df["AA"]==res, "#Codon"].tolist()
    freq_list = df.loc[df["AA"]==res, "Fraction"].tolist()
    codon_dict[res] = codon_list, freq_list


# ### Back-translate and codon optimize

target_dna = [] # List containing target DNA name/sequence pair(s)

for seq in target_msa:
    align_index = -1
    dna_index = 0
    seq_name = seq[0]
    seq_aa = seq[1]
    seq_dna = ""
    
    for res in seq_aa:
        align_index += 1        
        
        if template_aa[align_index] != "-":
            if res == "-": 
                dna_index += 3                
                continue
            elif res == template_aa[align_index]:
                if (re.search("TAA|TAG|TGA", seq_dna[-3:])) and (res != "M"): # Add start codon if missing
                    seq_dna += "ATG"                   
                codon = template_dna[dna_index:dna_index+3]
                seq_dna += codon
                dna_index += 3 
            else:
                if (re.search("TAA|TAG|TGA", seq_dna[-3:])) and (res != "M"): # Add start codon if missing
                    seq_dna += "ATG"
                rand_codon = random.choices(codon_dict[res][0], codon_dict[res][1], k=1)
                seq_dna += rand_codon[0]
                dna_index += 3 
        else:
            if res == "-":
                continue
            else:
                if (re.search("TAA|TAG|TGA", seq_dna[-3:])) and (res != "M"): # Add start codon if missing
                    seq_dna += "ATG"              
                rand_codon = random.choices(codon_dict[res][0], codon_dict[res][1], k=1)
                seq_dna += rand_codon[0]
                
    target_dna.append(tuple([seq_name, seq_dna]))


# Write FASTA file with optimized DNA sequences

out_f = os.path.splitext(msa_f)[0]+".codon_opt.fasta"

with open(out_f, "wt") as newfile:
    for seq in target_dna:
        newfile.write(">"+seq[0]+"\n")
        newfile.write(seq[1]+"\n")
        
print("Done. " + out_f + " contains codon-optimized sequences.")

