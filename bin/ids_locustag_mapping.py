#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 10:41:17 2025

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

        best_hit_id = None
        length_seq = None
        total_aligned = 0  # suma de dominios
        evalue = None
        found_first_hit = False

        # Parse hmmsearch result file
        with open(result_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.split()
                this_hit = parts[0]

                # Caso 1: primera línea útil → este es el mejor hit
                if not found_first_hit:
                    best_hit_id = this_hit
                    found_first_hit = True
                    #length_seq = float(parts[2])      # largo de la secuencia
                    evalue = float(parts[6])

                # Caso 2: ya tengo el mejor hit, solo acumulo dominios del mismo hit
                if this_hit == best_hit_id:
                    ali_from = int(parts[17])
                    ali_to = int(parts[18])
                    domain_len = ali_to - ali_from + 1
                    total_aligned += int(domain_len)
                else:
                    # llegó otro hit → dejo de leer, ya no me importa nada más
                    break

        # Si no hubo hits
        if not found_first_hit:
            return ("no hits")
        
        #coverage = total_aligned / length_seq

        # search GBFF for matching protein id
        # hit_id_clean = best_hit_id if best_hit_id else None
        
        if not best_hit_id:
            return None

        records = list(SeqIO.parse(gbff_file, "genbank"))
        for record in records:
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                qualifiers = feature.qualifiers
                protein_id = qualifiers.get("protein_id", [""])[0]
                locus_tag = qualifiers.get("locus_tag", [""])[0]
                # all_ids = [protein_id.split(".")[0]]  # keep same logic as before
                all_ids = [protein_id]  # keep same logic as before
                # if hit_id_clean in all_ids:
                if best_hit_id in all_ids:
                    found_locus_tag = locus_tag if locus_tag else protein_id
                    sequence = qualifiers.get("translation", [""])[0]
                    length_seq = len(sequence)
                    # return ("mapeado", hit_id_clean, found_locus_tag, length_seq, total_aligned, evalue)
                    return ("mapeado", best_hit_id, found_locus_tag, length_seq, total_aligned, evalue)

        # not found in GBFF
        # return ("no_mapeado", hit_id_clean)
        return ("no_mapeado", best_hit_id)

    except Exception as e:
        print(f"Error processing {result_path}: {e}", file=sys.stderr)
        return None

if __name__ == "__main__":
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(no_mapeadas_file) or ".", exist_ok=True)
    resultado = procesar_archivo(result_file)
    if resultado:
        if resultado[0] == "mapeado":
            _, protein_id, locus_tag, length_seq, total_aligned, evalue = resultado
            cog = result_file.split("_")[1].split(".")[0] # example_COG.result
            with open(output_file, "a") as out:
                out.write(f"{protein_id} -> {locus_tag}\tcog={cog}\tlen={length_seq}\taligned={total_aligned}\tevalue={evalue}\n")
        elif resultado == "no hits":
            with open(output_file, "a") as out:
                out.write(f"no hits\n")
        elif resultado[0] == "no_mapeado":
            # _, hit_id_clean = resultado
            _, best_hit_id = resultado
            with open(no_mapeadas_file, "a") as out:
                # out.write(f"{hit_id_clean}\n")
                out.write(f"{best_hit_id}\n")
    else:
        # nothing to write (no hits or errors logged to stderr)
        pass