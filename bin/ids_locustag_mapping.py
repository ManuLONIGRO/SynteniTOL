#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 10:41:17 2025

Create dataframe of synteny fwda - Wood Ljundahl pathway

UPDATED TO NEXTFLOW

Headers format:
taxid|species_name|isolate|n_repetitions

@author: mlonigro
"""


# Script para mapear IDs de proteínas encontradas por hmmsearch a sus locus_tag en archivos .gbff.
# Procesa múltiples archivos .result en paralelo para acelerar el análisis.
# Guarda los mapeos en 'id_locustag.txt' y los IDs no encontrados en 'ids_no_mapeadas.txt'.

import os
from Bio import SeqIO
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--genomes_dir", required=True, help="Path to the directory containing .gbff files")
parser.add_argument("--result_file", required=True, help="Path to a single .result file")
parser.add_argument("--output_file", required=True, help="Path to the output mapping file (append mode)")
parser.add_argument("--no_mapped_file", required=True, help="Path to the output file for unmapped IDs (append mode)")
args = parser.parse_args()

genomes_dir = args.genomes_dir
result_file = args.result_file
output_file = args.output_file
no_mapeadas_file = args.no_mapped_file

# Parameters
evalue_limit = 1e-5

def procesar_archivo(result_path):
    try:
        if not os.path.exists(result_path):
            print(f"Result file not found: {result_path}", file=sys.stderr)
            return None

        filename = os.path.basename(result_path)
        # original code assumed filename like "<organism>_... .result"
        organism_name = filename.split("_")[0]
        gbff_file = os.path.join(genomes_dir, f"{organism_name}.gbff")
        # name_file_gbff = f"{organism_name}.gbff"

        if not os.path.exists(gbff_file):
            print(f"GBFF not found for {organism_name}: {gbff_file}", file=sys.stderr)
            return None

        with open(result_path, "r") as f:
            lines = f.readlines()

        # quick no-hits check (keep original index logic)
        if len(lines) >= 18 and "No hits detected" in lines[15]:
            return None

        hit_presence = False
        hit_id = None
        evalue = None

        if len(lines) >= 15:
            columns = lines[14].split()
            if len(columns) > 8:
                try:
                    evalue = float(columns[0])
                    if evalue < evalue_limit:
                        hit_id = columns[8].split("|")[-1]
                        hit_presence = True
                except Exception:
                    pass

        if not hit_presence:
            return None

        # compute total aligned length from domain annotation (as in original)
        total_aligned_length = 0
        for i, line in enumerate(lines):
            if "Domain annotation for each sequence (and alignments):" in line:
                line_hits = i + 4
                break
        else:
            # no domain annotation section -> try to continue with default 0
            line_hits = None

        if line_hits is not None:
            max_domains_to_check = 10
            for offset in range(max_domains_to_check):
                idx = line_hits + offset
                if idx >= len(lines):
                    break
                fields = lines[idx].strip().split()
                if not fields or not fields[0].isdigit():
                    break
                try:
                    ali_from = int(fields[9])
                    ali_to = int(fields[10])
                    aligned_length = ali_to - ali_from + 1
                    total_aligned_length += aligned_length
                except Exception:
                    # ignore malformed domain lines
                    continue

        # search GBFF for matching protein id
        hit_id_clean = hit_id.split(".")[0] if hit_id else None
        if not hit_id_clean:
            return None

        records = list(SeqIO.parse(gbff_file, "genbank"))
        for record in records:
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                qualifiers = feature.qualifiers
                protein_id = qualifiers.get("protein_id", [""])[0]
                locus_tag = qualifiers.get("locus_tag", [""])[0]
                all_ids = [protein_id.split(".")[0]]  # keep same logic as before
                if hit_id_clean in all_ids:
                    found_id = locus_tag if locus_tag else protein_id
                    sequence = qualifiers.get("translation", [""])[0]
                    length_seq = len(sequence)
                    return ("mapeado", hit_id_clean, found_id, length_seq, total_aligned_length, evalue)

        # not found in GBFF
        return ("no_mapeado", hit_id_clean)

    except Exception as e:
        print(f"Error processing {result_path}: {e}", file=sys.stderr)
        return None

if __name__ == "__main__":
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(no_mapeadas_file) or ".", exist_ok=True)

    resultado = procesar_archivo(result_file)
    if resultado:
        if resultado[0] == "mapeado":
            _, seq_id, final_id, length_seq, total_aligned_length, evalue = resultado
            with open(output_file, "a") as out:
                out.write(f"{seq_id} -> {final_id}\tlen={length_seq}\taligned={total_aligned_length}\tevalue={evalue}\n")
        elif resultado[0] == "no_mapeado":
            _, hit_id_clean = resultado
            with open(no_mapeadas_file, "a") as out:
                out.write(f"{hit_id_clean}\n")
    else:
        # nothing to write (no hits or errors logged to stderr)
        pass