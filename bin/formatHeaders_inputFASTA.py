#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Nov 12 09:23 2025

Take the inputFASTA of NEXTFLOW and format the headers using dataframe.csv

@author: mlonigro
"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input_fasta", type=str, required=True, help="Input FASTA file path")
parser.add_argument("--dataframe", type=str, required=True, help="Path to the dataframe CSV file")
parser.add_argument("--output_fasta", type=str, required=True, help="Output formatted FASTA file path")
args = parser.parse_args()

input_fasta = args.input_fasta
dataframe = args.dataframe
output_fasta = args.output_fasta

# Load dataframe
df = pd.read_csv(dataframe, sep="\t")

# Build a mapping from protein_id to (taxid, species_name, isolate)
protein_id_to_info = {}
for _, row in df.iterrows():
    taxid = row["taxid"]
    species_name = row["specie"].replace(" ", "_")
    isolate = row["isolate"]
    protein_id = row["protein_id"]
    protein_id_to_info[protein_id] = (taxid, species_name, isolate)

# create header to replace for protein_id, create n_repetition of the header to add |n_repetitions to the new_header
header_count = {}
with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            protein_id = line.split('.')[0].strip().replace(">","")
            print(protein_id)
            if protein_id in protein_id_to_info:
                taxid, species_name, isolate = protein_id_to_info[protein_id]
                # Count repetitions
                if protein_id not in header_count:
                    header_count[protein_id] = 0
                header_count[protein_id] += 1
                n_repetitions = header_count[protein_id]
                new_header = f">{taxid}|{species_name}|{isolate}|{n_repetitions}\n"
                outfile.write(new_header)
            else:
                # If protein_id not found, keep original header
                outfile.write(line)
        else:
            outfile.write(line)

