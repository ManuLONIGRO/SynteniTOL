#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 13:36:00 2025

UPDATED FOR NEXTFLOW

@author: mlonigro
"""
import os
import csv
import shutil
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--tsv", required=True, help="Path to input TSV file with assembly info")
parser.add_argument("--result_zip_ncbi_dir", required=True, help="Path to the result of .zip ncbi_dataset download")
parser.add_argument("--output_dir_genomes", required=True, help="Path to output directory for renamed .gbff files")
args = parser.parse_args()

tsv_path = args.tsv
base_dir = args.result_zip_ncbi_dir
output_dir = args.output_dir_genomes


#tsv_path = "/home/mlonigro/pvc/blast/test/protein_to_organism_map.tsv"
#base_dir = "/home/mlonigro/pvc/blast/test/small_gbffs/ncbi_dataset/data/"       # donde están las carpetas de cada assembly
#output_dir = "/home/mlonigro/pvc/blast/test/genomes"    # donde querés guardar los .gbff copiados

os.makedirs(output_dir, exist_ok=True)

with open(tsv_path, newline="") as tsvfile:
    reader = csv.DictReader(tsvfile, delimiter="\t")
    for row in reader:
        assembly = row["assembly"].strip()
        protein_id = row["protein"].strip()
        # organism = row["organism"].strip().lower().replace(" ", "").replace("-","").replace(":","")
        organism = re.sub(r"[^a-zA-Z0-9]", "", row["organism"].strip().lower())
        # isolate = row["isolate"].strip().lower().replace(" ", "") if "isolate" in row and row["isolate"].strip() else ""
        raw_isolate = row.get("isolate", "")
        isolate = re.sub(r"[^a-zA-Z0-9]", "", raw_isolate.strip().lower()) if raw_isolate.strip() else ""

        src_dir = os.path.join(base_dir, assembly)
        if not os.path.isdir(src_dir):
            print(f"[!] The dir: {src_dir} it doesn't exists.")
            continue

        # buscar archivo .gbff dentro del directorio
        gbff_files = [f for f in os.listdir(src_dir) if f.endswith(".gbff")]
        if not gbff_files:
            print(f"[!] There isn't .gbff in {src_dir}")
            continue

        src_file = os.path.join(src_dir, gbff_files[0])
        if isolate:
            new_name = f"{organism}{isolate}{protein_id}.gbff".replace("_","")
        else:
            new_name = f"{organism}{protein_id}.gbff".replace("_","")

        dst_file = os.path.join(output_dir, new_name)
        shutil.copy2(src_file, dst_file)
        print(f"[+] Copy {src_file} -> {dst_file}")
