#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Nov 12 09:23 2025

Take the inputFASTA of NEXTFLOW and format the headers using dataframe.csv

@author: mlonigro
"""

import argparse
import pandas as pd
import os
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--input_fasta", type=str, required=True, help="Input FASTA file path")
#parser.add_argument("--dataframe", type=str, required=True, help="Path to the dataframe CSV file")
parser.add_argument("--output_fasta", type=str, required=True, help="Output formatted FASTA file path")
parser.add_argument("--genomes_dir", type=str, required=False, help="Path to the directory containing .gbff files")
# parser.add_argument("--protein_to_organism_map_tsv", type=str, required=True, help="Path to the protein to organism mapping TSV file")
args = parser.parse_args()

input_fasta = args.input_fasta
#dataframe = args.dataframe
# protein_to_organism_map_tsv = args.protein_to_organism_map_tsv
output_fasta = args.output_fasta
genomes_dir = args.genomes_dir

# Get the protein_id from the fasta file. Format is >protein_id:other info
#protein_ids = set()
protein_ids = []
with open(input_fasta, "r") as infile:
    for line in infile:
        if line.startswith(">"):
            protein_id = line.split(':')[0].strip().replace(">","")
            protein_ids.append(protein_id)

# Read gbff looking for taxid, species_name, isolate for each protein_id in protein_ids list
protein_id_to_info = {}
for gbff_file in os.listdir(genomes_dir):
    if gbff_file.endswith(".gbff"):
        gbff_path = os.path.join(genomes_dir, gbff_file)
        for record in SeqIO.parse(gbff_path, "genbank"):
            # source = record.annotations.get("source","")
            # organism = record.annotations.get("organism","")
            organism = "NA"
            isolate = "NA"
            strain = "NA"
            taxid = "NA"
            for feature in record.features:
                if feature.type == "source":
                    db_xrefs = feature.qualifiers.get("db_xref", [])
                    for xref in db_xrefs:
                        if xref.startswith("taxon:"):
                            taxid = xref.split(":")[1]
                            break
                    organism = feature.qualifiers.get("organism", ["NA"])[0].replace(" ", "_")
                    strain = feature.qualifiers.get("strain", ["NA"])[0].replace(" ", "_").replace(";","")
                    isolate = feature.qualifiers.get("isolate", ["NA"])[0].replace(" ", "_").replace(";","")

                    # # taxid
                    # if "db_xref" in feature.qualifiers:
                    #     for x in feature.qualifiers["db_xref"]:
                    #         if x.startswith("taxon:"):
                    #             taxid = x.split(":")[1]
                    # # isolate
                    # if "isolate" in feature.qualifiers:
                    #     isolate = feature.qualifiers["isolate"][0]
                    
                # Search protein_id in CDS
                if feature.type == "CDS" and "protein_id" in feature.qualifiers:
                    pid = feature.qualifiers["protein_id"][0]
                    if pid in protein_ids:
                        protein_id_to_info[pid] = {
                            "taxid": (taxid if taxid is not None else "NA"),
                            "organism": (organism if organism is not None else "NA"),
                            "strain": (strain if strain is not None else "NA"),
                            "isolate": (isolate if isolate is not None else "NA"),                            
                        }

# create header to replace for protein_id, create n_repetition of the header to add |n_repetitions to the new_header
header_count = {}
with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            protein_id = line.split(':')[0].strip().replace(">","")
            if protein_id in protein_id_to_info:
                info = protein_id_to_info[protein_id]
                taxid = info.get("taxid", "")
                species_name = info.get("organism", "NA").replace(" ", "_")
                strain = info.get("strain", "NA").replace(" ", "_").replace(";","")
                isolate = info.get("isolate", "NA").replace(" ", "_").replace(";","")
                # Count repetitions
                header = f">{taxid}|{species_name}|{strain}|{isolate}"
                if header not in header_count:
                    header_count[header] = 0
                header_count[header] += 1
                n_repetitions = header_count[header]
                new_header = f"{header}|{n_repetitions}\n"
                outfile.write(new_header)
            else:
                # If protein_id not found, keep original header
                outfile.write(line)
        else:
            outfile.write(line)

