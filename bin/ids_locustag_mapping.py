#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 10:41:17 2025

UPDATED TO NEXTFLOW

Headers format:
taxid|species_name|isolate|n_repetitions

@author: mlonigro
"""
# Script to map protein IDs found by hmmsearch to their locus_tag in .gbff files.
# Processes multiple .result files in parallel to speed up the analysis.
# Saves the mappings in 'id_locustag.txt' and the IDs not found in 'ids_no_mapeadas.txt'.

import os
from Bio import SeqIO
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--genomes_dir", required=True, help="Path to the directory containing .gbff files")
parser.add_argument("--result_file", required=True, help="Path to a single .result file")
parser.add_argument("--gene_name", required=False, help="Gene/profile label to annotate mapping rows")
parser.add_argument("--output_file", required=True, help="Path to the output mapping file (append mode)")
parser.add_argument("--no_mapped_file", required=True, help="Path to the output file for unmapped IDs (append mode)")
parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold for filtering")
args = parser.parse_args()

genomes_dir = args.genomes_dir
result_file = args.result_file
output_file = args.output_file
no_mapeadas_file = args.no_mapped_file
evalue_limit = args.evalue

def procesar_archivo(result_path):
    try:
        if not os.path.exists(result_path):
            print(f"Result file not found: {result_path}", file=sys.stderr)
            return []

        filename = os.path.basename(result_path)
        organism_name = filename.split("_")[0]
        gbff_file = os.path.join(genomes_dir, f"{organism_name}.gbff")

        if not os.path.exists(gbff_file):
            print(f"GBFF not found for {organism_name}: {gbff_file}", file=sys.stderr)
            return None

        # Dictionary to accumulate information for each hit with e‑value <= limit
        hit_stats = {}

        # Parse hmmsearch result file and save ALL hits with e‑value <= evalue_limit
        with open(result_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.split()
                if len(parts) < 19:
                    # strange/incomplete line; skip it
                    continue

                # get all hits with an evalue <= evalue_limit
                evalue = float(parts[6])
                if evalue > evalue_limit:
                    continue

                hit_id = parts[0]
                ali_from = int(parts[17])
                ali_to = int(parts[18])
                domain_len = ali_to - ali_from + 1
                stats = hit_stats.setdefault(hit_id, {"total_aligned": 0, "evalue": evalue})
                stats["total_aligned"] += int(domain_len)
                
                # Get the acc
                hmm_acc = float(parts[21])
                print(f"The accuracy for hmmsearch of {hit_id} is {hmm_acc}")
                stats["hmm_acc"] = hmm_acc
                
                # Save the best (minimum) e‑value observed for that hit
                if evalue < stats["evalue"]:
                    stats["evalue"] = evalue
                
                

        # If there were no hits below the limit
        if not hit_stats:
            return [("no hits",)]

        resultados = []

        records = list(SeqIO.parse(gbff_file, "genbank"))
        for record in records:
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                qualifiers = feature.qualifiers
                protein_id = qualifiers.get("protein_id", [""])[0]
                if not protein_id:
                    continue

                # We keep the ID as it appears in hmmsearch (without the suffix)
                if protein_id not in hit_stats:
                    continue

                locus_tag = qualifiers.get("locus_tag", [""])[0]
                found_locus_tag = locus_tag if locus_tag else protein_id
                sequence = qualifiers.get("translation", [""])[0]
                length_seq = len(sequence)

                stats = hit_stats[protein_id]
                total_aligned = stats["total_aligned"]
                evalue = stats["evalue"]
                hmmsearch_acc = stats["hmm_acc"]

                resultados.append(
                    ("mapeado", protein_id, found_locus_tag, length_seq, total_aligned, evalue, hmmsearch_acc)
                )

        # Any hit that remained in hit_stats but was not found in the GBFF is marked as not mapped
        mapped_ids = {r[1] for r in resultados if r[0] == "mapeado"}
        for hit_id in hit_stats.keys():
            if hit_id not in mapped_ids:
                resultados.append(("no_mapeado", hit_id))

        return resultados

    except Exception as e:
        print(f"Error processing {result_path}: {e}", file=sys.stderr)
        return []

if __name__ == "__main__":
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(no_mapeadas_file) or ".", exist_ok=True)
    resultado = procesar_archivo(result_file)
    if not resultado:
        # Nothing to write (no hits or errors logged to sys.stderr)
        sys.exit(0)

    # Prefer explicit label from workflow; fallback keeps backward compatibility
    cog = args.gene_name
    if not cog:
        try:
            cog = result_file.split("_")[1].split(".")[0]
        except Exception:
            cog = "unknown"

    with open(output_file, "a") as out_map, open(no_mapeadas_file, "a") as out_no_map:
        for item in resultado:
            etiqueta = item[0]
            if etiqueta == "mapeado":
                _, protein_id, locus_tag, length_seq, total_aligned, evalue, hmmsearch_acc  = item
                out_map.write(
                    f"{protein_id} -> {locus_tag}\tcog={cog}\tlen={length_seq}\t"
                    f"aligned={total_aligned}\tevalue={evalue}\thmmsearch_acc={hmmsearch_acc}\n"
                )
            elif etiqueta == "no hits":
                out_map.write("no hits\n")
            elif etiqueta == "no_mapeado":
                _, hit_id = item
                out_no_map.write(f"{hit_id}\n")