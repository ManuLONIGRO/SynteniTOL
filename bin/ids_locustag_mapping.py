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
        # original code assumed filename like "<organism>_... .result"
        organism_name = filename.split("_")[0]
        gbff_file = os.path.join(genomes_dir, f"{organism_name}.gbff")
        # name_file_gbff = f"{organism_name}.gbff"

        if not os.path.exists(gbff_file):
            print(f"GBFF not found for {organism_name}: {gbff_file}", file=sys.stderr)
            return None

        # Diccionario para acumular información por cada hit con e‑value <= límite
        # hit_stats[hit_id] = {"total_aligned": int, "evalue": float}
        hit_stats = {}

        # Parse hmmsearch result file y guardar TODOS los hits con e‑value <= evalue_limit
        with open(result_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.split()
                if len(parts) < 19:
                    # línea rara / incompleta; la saltamos
                    continue

                # obtener todos los hits que tengan un evalue <= evalue_limit
                evalue = float(parts[6])
                if evalue > evalue_limit:
                    continue

                hit_id = parts[0]
                ali_from = int(parts[17])
                ali_to = int(parts[18])
                domain_len = ali_to - ali_from + 1
                stats = hit_stats.setdefault(hit_id, {"total_aligned": 0, "evalue": evalue})
                stats["total_aligned"] += int(domain_len)
                
                # Get the coverage acc
                hmm_coverage = float(parts[21])
                print(f"The coverage for hmmsearch of {hit_id} is {hmm_coverage}")
                stats["hmm_coverage"] = hmm_coverage
                
                # Guardamos el mejor (mínimo) e‑value observado para ese hit
                if evalue < stats["evalue"]:
                    stats["evalue"] = evalue
                
                

        # Si no hubo hits por debajo del límite
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

                # Nos quedamos con el ID tal cual aparece en hmmsearch (sin recortar el sufijo)
                if protein_id not in hit_stats:
                    continue

                locus_tag = qualifiers.get("locus_tag", [""])[0]
                found_locus_tag = locus_tag if locus_tag else protein_id
                sequence = qualifiers.get("translation", [""])[0]
                length_seq = len(sequence)

                stats = hit_stats[protein_id]
                total_aligned = stats["total_aligned"]
                evalue = stats["evalue"]
                hmmsearch_cover = stats["hmm_coverage"]

                resultados.append(
                    ("mapeado", protein_id, found_locus_tag, length_seq, total_aligned, evalue, hmmsearch_cover)
                )

        # Cualquier hit que quedó en hit_stats pero no se encontró en el GBFF se marca como no mapeado
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
        # nothing to write (no hits or errors logged to stderr)
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
                _, protein_id, locus_tag, length_seq, total_aligned, evalue, hmmsearch_cover  = item
                out_map.write(
                    f"{protein_id} -> {locus_tag}\tcog={cog}\tlen={length_seq}\t"
                    f"aligned={total_aligned}\tevalue={evalue}\thmmsearch_cover={hmmsearch_cover}\n"
                )
            elif etiqueta == "no hits":
                out_map.write("no hits\n")
            elif etiqueta == "no_mapeado":
                _, hit_id = item
                out_no_map.write(f"{hit_id}\n")