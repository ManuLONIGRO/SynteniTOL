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
import time
import os

parser = argparse.ArgumentParser()
parser.add_argument("--fasta",required=True, help="Path to input FASTA file")
parser.add_argument("--out_tsv",required=True, help="Path to  output TSV file")
parser.add_argument("--out_assemblies", required=True, help="Path to output assemblies file")
#parser.add_argument("--out_taxonomy", required=True, help="Path to output taxonomy TSV file")
parser.add_argument("--out_no_assembly_list", required=True, help="Path to output nucleotide accessions with no assembly info")
# parser.add_argument("--ncbi_api_key", required=False, help="NCBI API key to increase request limits")
args = parser.parse_args()

fasta_path = args.fasta
output_tsv = args.out_tsv
output_assemblies = args.out_assemblies
#output_taxonomy_tsv = args.out_taxonomy
output_no_assembly_list = args.out_no_assembly_list

# if args.ncbi_api_key:
#     os.environ["NCBI_API_KEY"] = args.ncbi_api_key

protein_id_list = []
with open(fasta_path) as fasta:
    for line in fasta:
        if line.startswith(">"):
            protein_id = line.split(":")[0].replace(">", "").strip()
            protein_id_list.append(protein_id)

# Convert list to format suitable for efetch
queries = ",".join(protein_id_list)
print("START PARSING", flush=True)
# Execute efetch (uses ipg format)
cmd_ipg = f"efetch -db protein -id {queries} -format ipg"
result = subprocess.check_output(cmd_ipg, shell=True, text=True)

# Process lines and save in list
rows = []
for line in result.splitlines():
    if not line or line.startswith("Id"):
        continue  # skip header or empty lines
    parts = line.split("\t")
    # if len(parts) == 11:
    while len(parts) < 11:
        parts.append("NA")
    
    row = {
        #ipg_id=parts[0], COULD BE INTERESTING, ID FOR THE SEQUENCE, CAN BE PRESENT IN SEVERAL EQUAL SEQUENCES OR ENTRIES
        #source=parts[1],
        "nucleotide_accession": parts[2] or "NA",
        #start=parts[3],
        #stop=parts[4],
        #strand=parts[5],
        "protein": parts[6] or "NA",
        #protein_name=parts[7],
        "organism": parts[8] or "NA",
        "strain": parts[9].replace(";","") or "NA",
        "assembly": parts[10] or "NA",
    }
    
    rows.append(row)

# When searching in ipg, I get all entries that have the same sequence as my protein_id, filter by the protein_id to get only the entries that are in the fasta, and not download more gbff.
rows = [r for r in rows if r["protein"] in protein_id_list]

# Eliminate duplicates based on protein ID
seen = set()
new_rows = []

for r in rows:
    pid = r["protein"]
    if pid not in seen:
        seen.add(pid)
        new_rows.append(r)
print("MIDDLE PARSING", flush=True)
# Save the assemblies for datasets --inputfile and nucleotide accessions with no assembly info
no_assembly_list = [r["nucleotide_accession"] for r in new_rows if r["assembly"] in ("NA", None, "")]
assemblies = {r["assembly"] for r in new_rows if r["assembly"] and r["assembly"]!= "NA"}
# assemblies = {r["assembly"] for r in new_rows if r["assembly"]} # Create a set to remove duplicates

# Save assemblies to download in a file
with open(output_assemblies,"w") as file_assemblies:
    for a in assemblies:
        file_assemblies.write(a + "\n")

# Save nucleotide accessions with no assembly info
with open(output_no_assembly_list, "w") as no_assembly_file:
    for item in no_assembly_list:
        no_assembly_file.write(item + "\n")
print("END PARSING", flush=True)
# Execute efetch only if we have assemblies to query
assembly_to_isolate = {}
if assemblies:
    assemblies_query = ",".join(assemblies)
    cmd_assembly = f"efetch -db assembly -id {assemblies_query} -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession Sub_value Isolate"
    result2 = subprocess.check_output(cmd_assembly, shell=True, text=True, timeout=20).strip()

    # Diccionario assembly: isolate
    for line2 in result2.splitlines():
        parts2 = line2.strip().split("\t")
        acc = parts2[0] if len(parts2) > 0 else "NA"
        isolate = parts2[1].replace(" ", "_") if len(parts2) > 1 else "NA"
        assembly_to_isolate[acc] = isolate

# Add the isolate to each row
for r in new_rows:
    r["isolate"] = assembly_to_isolate.get(r["assembly"], "NA")

# # Clean protein_ids in rows
# protein_ids_extracted = []
# for r in rows:
#     pid_row = r["protein"]
#     if pid_row not in protein_ids_extracted:
#         protein_ids_extracted.append(pid_row)
#     else:
#         continue
# new_rows = []
# for pid in protein_ids_extracted:
#     for r in rows:
#         if r["protein"] == pid:
#             new_rows.append(r)
#             break  # only the first occurrence

fieldnames = ["protein", "organism", "strain", "assembly", "isolate", "nucleotide_accession"]

# Save in TSV with header
with open(output_tsv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(new_rows)
#----------------------------------------------------------------------------
# # For each assembly, get taxonomy info
# with open(output_taxonomy_tsv, "a") as tax_file:
#             tax_file.write(f"assembly\ttaxid\tscientific_name\trank\n") # The first line of the file
# if assemblies:
#     for accession in assemblies:
#         cmd_taxonomy = f"esearch -db assembly -query {accession} \
#             | elink -target taxonomy \
#             | efetch -format xml \
#             | xtract -pattern LineageEx -group Taxon -sep '\t' -tab '\n' -element TaxId,ScientificName,Rank \
#             | awk -v acc='{accession}' 'BEGIN{{OFS=\"\t\"}} {{print acc, $0}}'"
#         try:
#             result_taxonomy = subprocess.check_output(cmd_taxonomy, shell=True, text=True, timeout= 40).strip()
#         except subprocess.CalledProcessError as e:
#             result = ""
        
#         # Save taxonomy info in a file
#         for line in result_taxonomy.splitlines():
#             parts = line.strip().split("\t")
#             acc = parts[0] if len(parts) > 0 else "NA"
#             taxid = parts[1] if len(parts) > 1 else "NA"
#             scientific_name = parts[2].replace(" ","_") if len(parts) > 2 else "NA"
#             rank = parts[3].replace(" ","_") if len(parts) > 3 else "NA"
#             with open(output_taxonomy_tsv, "a") as tax_file:
#                 tax_file.write(f"{acc}\t{taxid}\t{scientific_name}\t{rank}\n")
#         time.sleep(0.34)  # To avoid overloading NCBI servers

