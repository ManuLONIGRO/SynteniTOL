#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Nextflow-friendly script to build the WLP dataframe TSV.

Inputs:
  --genomes_dir    Path to directory containing .gbff genomes
  --mappings_file  Path to merged mappings file (seq_id -> locus_tag\tlen=..\taligned=..\tevalue=..)
  --results_dir    Directory containing *.result files (to infer gene names)
  --protein_map    TSV protein->organism mapping (columns: protein, organism, strain, assembly, isolate)
  --out_tsv        Output TSV path

It reconstructs the dictionaries used by the original workflow and writes a
dataframe with columns:
  taxid, especie, isolate, locus_tag, id_proteina, gene, evalue, longitud, alineado, coverage
"""

import argparse
import os
from typing import Dict, Tuple, List
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes_dir", required=True)
    parser.add_argument("--mappings_file", required=True)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--protein_map", required=True)
    parser.add_argument("--out_tsv", required=True)
    return parser.parse_args()


def load_id_to_gene_from_results(results_dir: str) -> Dict[str, str]:
    id_to_gene: Dict[str, str] = {}
    for name in os.listdir(results_dir):
        if not name.endswith(".result"):
            continue
        # Gene name from filename, last token before .result
        parts = name.split("_")
        if len(parts) < 2:
            continue
        gene_name = parts[-1].replace(".result", "")

        # Extract best-hit protein id from line 15 (index 14)
        result_path = os.path.join(results_dir, name)
        try:
            with open(result_path, "r") as fh:
                lines = fh.readlines()
            if len(lines) >= 18 and "No hits detected" in lines[15]:
                continue
            if len(lines) >= 15:
                cols = lines[14].split()
                if len(cols) > 8:
                    hit_id = cols[8].split("|")[-1].split(".")[0]
                    id_to_gene[hit_id] = gene_name
        except Exception:
            continue
    return id_to_gene


def load_mappings(mappings_file: str) -> Tuple[Dict[str, str], Dict[str, int], Dict[str, int], Dict[str, str]]:
    # Returns: id->locus, locus->len, locus->aligned, locus->evalue
    id_to_locus: Dict[str, str] = {}
    locus_to_length: Dict[str, int] = {}
    locus_to_aligned: Dict[str, int] = {}
    locus_to_evalue: Dict[str, str] = {}

    with open(mappings_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            # Format: <seq_id> -> <locus>\tlen=<len>\taligned=<n>\tevalue=<e>
            try:
                left, right = line.split(" -> ", 1)
                seq_id = left
                fields = right.split("\t")
                locus_tag = fields[0]
                id_to_locus[seq_id] = locus_tag
                # Parse metrics if present
                for field in fields[1:]:
                    if field.startswith("len="):
                        try:
                            locus_to_length[locus_tag] = int(field.split("=", 1)[1])
                        except Exception:
                            pass
                    elif field.startswith("aligned="):
                        try:
                            locus_to_aligned[locus_tag] = int(field.split("=", 1)[1])
                        except Exception:
                            pass
                    elif field.startswith("evalue="):
                        locus_to_evalue[locus_tag] = field.split("=", 1)[1]
            except ValueError:
                # line not conforming; skip
                continue

    return id_to_locus, locus_to_length, locus_to_aligned, locus_to_evalue


def load_protein_map(tsv_path: str) -> Dict[str, Tuple[str, str]]:
    mapping: Dict[str, Tuple[str, str]] = {}
    with open(tsv_path) as fh:
        header = fh.readline().strip().split("\t")
        # expected: protein, organism, strain, assembly, isolate
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or len(parts) < 2:
                continue
            protein_id = parts[idx.get("protein", 0)]
            organism = parts[idx.get("organism", 1)].replace(" ", "_")
            isolate = parts[idx.get("isolate", 4)] if len(parts) > 4 else "NA"
            mapping[protein_id.split(".")[0]] = (organism, isolate)
    return mapping


def iterate_genomes(genomes_dir: str):
    for name in os.listdir(genomes_dir):
        if not name.endswith(".gbff"):
            continue
        yield os.path.join(genomes_dir, name)


def build_dataframe(genomes_dir: str,
                    id_to_locus: Dict[str, str],
                    id_to_gene: Dict[str, str],
                    locus_to_length: Dict[str, int],
                    locus_to_aligned: Dict[str, int],
                    locus_to_evalue: Dict[str, str],
                    protein_map: Dict[str, Tuple[str, str]],
                    out_tsv: str) -> None:
    locus_to_id = {v: k for k, v in id_to_locus.items()}

    with open(out_tsv, "w") as out:
        out.write("taxid\tespecie\tisolate\tlocus_tag\tid_proteina\tgene\tevalue\tlongitud\talineado\tcoverage\n")

        for gbff_path in iterate_genomes(genomes_dir):
            records = list(SeqIO.parse(gbff_path, "genbank"))
            if not records:
                continue
            taxid = None
            for record in records:
                for feature in record.features:
                    if feature.type == "source":
                        db_xrefs = feature.qualifiers.get("db_xref", [])
                        for xref in db_xrefs:
                            if xref.startswith("taxon:"):
                                taxid = xref.split(":")[1]
                                break
                        if taxid:
                            break
                if taxid:
                    break
            if not taxid:
                continue

            # Emit rows for features whose locus_tag is in our mapping
            for record in records:
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    if not locus_tag:
                        continue
                    if locus_tag not in locus_to_id:
                        continue
                    protein_id = locus_to_id.get(locus_tag, "NA")
                    gene = id_to_gene.get(protein_id, "NA")
                    especie, isolate = protein_map.get(protein_id, ("Unknown", "NA"))
                    longitud = int(locus_to_length.get(locus_tag, 0))
                    aligned = int(locus_to_aligned.get(locus_tag, 0))
                    evalue = locus_to_evalue.get(locus_tag, "NA")
                    coverage = round(aligned / longitud, 2) if longitud > 0 else "NA"

                    out.write(
                        f"{taxid}\t{especie}\t{isolate}\t{locus_tag}\t{protein_id}\t{gene}\t{evalue}\t{longitud}\t{aligned}\t{coverage}\n"
                    )


def main():
    args = parse_args()
    id_to_locus, locus_to_length, locus_to_aligned, locus_to_evalue = load_mappings(args.mappings_file)
    id_to_gene = load_id_to_gene_from_results(args.results_dir)
    protein_map = load_protein_map(args.protein_map)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    build_dataframe(
        genomes_dir=args.genomes_dir,
        id_to_locus=id_to_locus,
        id_to_gene=id_to_gene,
        locus_to_length=locus_to_length,
        locus_to_aligned=locus_to_aligned,
        locus_to_evalue=locus_to_evalue,
        protein_map=protein_map,
        out_tsv=args.out_tsv,
    )


if __name__ == "__main__":
    main()


