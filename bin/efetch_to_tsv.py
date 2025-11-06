#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 15:17:13 2025

ALL QUERIES IN A LIST AND USE EFETCH TO ACCESS NCBI ONLY ONCE
Get the accession of the gbff of interest to not download by taxid that downloads more.

UPDATED FOR NEXTFLOW 28/10/2025
@author: mlonigro
"""
import argparse
import subprocess
import csv

parser = argparse.ArgumentParser()
parser.add_argument("--fasta",required=True, help="Path to input FASTA file")
parser.add_argument("--out-tsv",required=True, help="Path to  output TSV file")
parser.add_argument("--out-assemblies", required=True, help="Path to output assemblies file")
args = parser.parse_args()

fasta_path = args.fasta
output_tsv = args.out_tsv
output_assemblies = args.out_assemblies

protein_id_list = []
with open(fasta_path) as fasta:
    for line in fasta:
        if line.startswith(">"):
            protein_id = line.split(":")[0].replace(">", "").strip()
            protein_id_list.append(protein_id)

# Convert list to format suitable for efetch
queries = ",".join(protein_id_list)

# Execute efetch (uses ipg format)
cmd_ipg = f"efetch -db protein -id {queries} -format ipg"
result = subprocess.check_output(cmd_ipg, shell=True, text=True).strip()

# Process lines and save in list
rows = []
for line in result.splitlines():
    if not line.strip() or line.startswith("Id"):
        continue  # skip header or empty lines
    parts = line.strip().split("\t")
    if len(parts) == 11:
        row = dict(
            #ipg_id=parts[0], COULD BE INTERESTING, ID FOR THE SEQUENCE, CAN BE PRESENT IN SEVERAL EQUAL SEQUENCES OR ENTRIES
            #source=parts[1],
            #nucleotide_accession=parts[2],
            #start=parts[3],
            #stop=parts[4],
            #strand=parts[5],
            protein=parts[6],
            #protein_name=parts[7],
            organism=parts[8],
            strain=parts[9],
            assembly=parts[10],
        )
        rows.append(row)

# When searching in ipg, I get all entries that have the same sequence as my protein_id, filter by the protein_id to get only the entries that are in the fasta, and not download more gbff.
rows = [r for r in rows if r["protein"] in protein_id_list]

# Save the assemblies for datasets --inputfile
assemblies = {r["assembly"] for r in rows if r["assembly"]} # Create a set to remove duplicates
with open(output_assemblies,"w") as file_assemblies:
    for a in assemblies:
        file_assemblies.write(a + "\n")

assemblies_query = ",".join(assemblies)
cmd_assembly = f"efetch -db assembly -id {assemblies_query} -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession Sub_value Isolate"
result2 = subprocess.check_output(cmd_assembly, shell=True, text=True).strip()

# Diccionario assemlby: isolate
assembly_to_isolate = {}
for line2 in result2.splitlines():
    parts2 = line2.strip().split("\t")
    acc = parts2[0] if len(parts2) > 0 else ""
    isolate = parts2[1] if len(parts2) > 1 else ""
    assembly_to_isolate[acc] = isolate

# Add the isolate to each row
for r in rows:
    r["isolate"] = assembly_to_isolate.get(r["assembly"], "")

fieldnames = ["protein", "organism", "strain", "assembly", "isolate"]

# Save in TSV with header
with open(output_tsv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

