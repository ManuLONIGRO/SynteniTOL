#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 28 09:33:15 2025

UPDATED TO NEXTFLOW

Creates the proteome file (.prot) of the .gbff

@author: mlonigro
"""

import os
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--genomes_dir", required=True, help="Path to the directory containing .gbff files")
parser.add_argument("--proteomes_dir", required=True, help="Path to the output directory for .prot files")
args = parser.parse_args()

input_dir = args.genomes_dir
output_dir = args.proteomes_dir

os.makedirs(output_dir, exist_ok=True)  # The dir exists or is created
def procesar_gbff(filename):
    input_file = os.path.join(input_dir, filename)
    output_file = os.path.join(output_dir, filename.replace(".gbff", ".prot"))

    # If it already exists, skip
    if os.path.exists(output_file):
        return f"⏭️  {filename} ya procesado."

    try:
        with open(output_file, "w") as out_f:
            for record in SeqIO.parse(input_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        qualifiers = feature.qualifiers
                        if "translation" in qualifiers:
                            protein_seq = qualifiers["translation"][0]
                            protein_id = qualifiers.get("protein_id", ["no_protein_id"])[0]
                            out_f.write(f">{protein_id}\n{protein_seq}\n")
        return f"✔️  {filename} procesado."
    except Exception as e:
        return f"❌  Error en {filename}: {e}"

if __name__ == "__main__":
    gbff_files = [f for f in os.listdir(input_dir) if f.endswith(".gbff")]

    with ProcessPoolExecutor() as executor:
        for resultado in executor.map(procesar_gbff, gbff_files):
            print(resultado)