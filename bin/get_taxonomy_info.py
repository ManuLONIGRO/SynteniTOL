#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 09:45

Open output_assemblies from efetch_to_tsv.py and get the taxonomic information for each assembly using efetch with the taxid. Save the output in a TSV file with the following columns: assembly_accession, taxid, superkingdom, phylum, class, order, family, genus, species. If no taxonomy information is found for a given taxid, save "NA" in the corresponding fields.

For Nextflow pipeline

@autor: mlonigro
"""

import argparse
import subprocess
import time
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input_all_accessions_list", required=True, help="Path to input assemblies file")
parser.add_argument("--output_taxonomy_tsv", required=True, help="Path to output taxonomy TSV file")
args = parser.parse_args()

print(f"START GETTING TAXONOMY INFO using {args.input_all_accessions_list}")

# The NCBI API key is injected into the process environment by Nextflow
# (NCBI_API_KEY env var / secret), never passed via the command line.
# Use the API key if present (shorter sleep) to avoid overloading NCBI servers.
user_ncbi_api_key = bool(os.environ.get("NCBI_API_KEY"))

# Create output taxonomy file
with open(args.output_taxonomy_tsv, "a") as tax_file:
            tax_file.write(f"assembly\ttaxid\tscientific_name\trank\n") # The first line of the file

# load txt with list of assemblies
with open(args.input_all_accessions_list) as accessions_list:
    for accession in accessions_list:
        accession = accession.strip()  # remove newline and surrounding whitespace
        if not accession:
            continue
        print(accession)
        cmd_taxonomy = f"esearch -db assembly -query {accession} \
            | elink -target taxonomy \
            | efetch -format xml \
            | xtract -pattern LineageEx -group Taxon -sep '\t' -tab '\n' -element TaxId,ScientificName,Rank \
            | awk -v acc='{accession}' 'BEGIN{{OFS=\"\t\"}} {{print acc, $0}}'"
        try:
            result_taxonomy = subprocess.check_output(cmd_taxonomy, shell=True, text=True).strip()
        except subprocess.CalledProcessError as e:
            print(f"Error fetching taxonomy for accession:{accession} - {e}")
            result = ""
            result_taxonomy = "NA\tNA\tNA\tNA"
        except subprocess.TimeoutExpired as e:
            print(f"Timeout expired for accession:{accession}")
            result = ""
            result_taxonomy = "NA\tNA\tNA\tNA"
        
        # Save taxonomy info in a file
        for line in result_taxonomy.splitlines():
            parts = line.strip().split("\t")
            acc = parts[0] if len(parts) > 0 else "NA"
            taxid = parts[1] if len(parts) > 1 else "NA"
            scientific_name = parts[2].replace(" ","_") if len(parts) > 2 else "NA"
            rank = parts[3].replace(" ","_") if len(parts) > 3 else "NA"
            with open(args.output_taxonomy_tsv, "a") as tax_file:
                tax_file.write(f"{acc}\t{taxid}\t{scientific_name}\t{rank}\n")
        if user_ncbi_api_key:
            print("Using your NCBI API KEY")
            time.sleep(0.15)
        else:
            time.sleep(0.35) # To avoid overloading NCBI servers