#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 10:52:59 2025

Creates the files for itol from dataframe.tsv

Refactored for Nextflow: no multiprocessing, no hardcoded paths. Uses CLI args.
"""

import os
import argparse
from collections import defaultdict, Counter
from typing import Dict, List, Tuple

import pandas as pd
from Bio import SeqIO
from ete3 import NCBITaxa

parser = argparse.ArgumentParser()
parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold for filtering")
parser.add_argument("--input_dir", type=str, required=True, help="Input directory containing .gbff files")
parser.add_argument("--dataframe", type=str, required=True, help="Path to the dataframe TSV file")
parser.add_argument("--out_prefix", type=str, required=True, help="Output prefix for generated files")
parser.add_argument("--gene_list", type=str, nargs='*', default=None, help="List of genes/COGs of interest (space-separated)")
parser.add_argument("--cogs", type=str, default=None, help="Comma-separated COG list, e.g., COG1152,COG1795")
args = parser.parse_args()

choose_evalue = args.evalue
formatted_evalue = f"{choose_evalue:.0e}"  # scientific notation for filenames

# Load dataframe and build lookups
print(f"Reading dataframe from {args.dataframe}")

df = pd.read_csv(args.dataframe, sep="\t")
# Filter by coverage
df_filtrado = df[df["coverage"] > 0.65]

# Build lookups
locus_to_gene: Dict[str, str] = dict(zip(df_filtrado["locus_tag"], df_filtrado["gene"]))
dic_taxid: Dict[str, str] = dict(zip(df_filtrado["locus_tag"], df_filtrado["taxid"]))

print(f"Creating presence_binary_data_fromTSV_{formatted_evalue}...")

ncbi = NCBITaxa()

# Determine COG of interest and its locus_tags before scanning GBFFs
def parse_cogs(cogs_str: str) -> List[str]:
    if not cogs_str:
        return []
    return [c.strip() for c in cogs_str.split(',') if c.strip()]

gene_list: List[str] = []
if args.gene_list:
    gene_list = list(args.gene_list)
if not gene_list and args.cogs:
    gene_list = parse_cogs(args.cogs)
if not gene_list:
    raise SystemExit("You must provide --gene_list or --cogs")

cog_of_interest = gene_list[0]
locus_tags_cog_list = df_filtrado[df_filtrado["gene"].str.lower() == cog_of_interest.lower()]["locus_tag"].unique().tolist()
locus_goi_set = set([lt.lower() for lt in locus_tags_cog_list])

def process_gbff(input_file: str, locus_to_gene: Dict[str, str], locus_goi_set: set):
    records = list(SeqIO.parse(input_file, "genbank"))
    if not records:
        return None

    best_contig_genes: List[Tuple[int, int, int, str, str]] = []
    best_contig_len = 0
    taxid = None
    isolate = "NA"

    # choose only the contig that contains the GOI; if multiple, pick the one with more genes
    for record in records:
        contig_genes: List[Tuple[int, int, int, str, str]] = []
        contains_goi = False
        local_taxid = None
        local_isolate = "NA"

        for feature in record.features:
            if feature.type == "source":
                db_xrefs = feature.qualifiers.get("db_xref", [])
                local_isolate = feature.qualifiers.get("isolate", ["NA"])[0].replace(" ", "_")
                for xref in db_xrefs:
                    if xref.startswith("taxon:"):
                        local_taxid = xref.split(":")[1]
                        break
            if feature.type == "CDS":
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                if not locus_tag:
                    continue
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = int(feature.location.strand)
                gene_name = locus_to_gene.get(locus_tag, locus_tag)
                contig_genes.append((start, end, strand, gene_name, locus_tag))
                if locus_tag.lower() in locus_goi_set:
                    contains_goi = True

        if contains_goi and len(contig_genes) > len(best_contig_genes):
            best_contig_genes = contig_genes
            best_contig_len = len(record.seq)
            taxid = local_taxid
            isolate = local_isolate

    if not best_contig_genes or not taxid:
        return None

    # Resolve species name from taxid
    try:
        lineage = ncbi.get_lineage(int(taxid))
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        species_taxid = [tid for tid in lineage if ranks[tid] == "species"]
        if species_taxid:
            nombre_org = names[species_taxid[0]]
            species_name = nombre_org.replace(" ", "_").replace("(", "").replace(")", "")
        else:
            species_name = "Unknown"
    except Exception:
        species_name = "Unknown"

    header = f"{species_taxid[0]}|{species_name}|{isolate}"
    return (header, best_contig_len, best_contig_genes)

# Sequential processing (Nextflow manages parallelism across tasks if needed)
files = [os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir) if f.endswith(".gbff")]

genomes_data_binary_raw = []
for gbff in files:
    result = process_gbff(gbff, locus_to_gene, locus_goi_set)
    if result is not None:
        genomes_data_binary_raw.append(result)

# Assign repetition numbers per header
header_count = defaultdict(int)
genomes_data_numbered = []
for header, genome_length, all_genes in genomes_data_binary_raw:
    header_count[header] += 1
    header_with_rep = f"{header}|{header_count[header]}"
    genomes_data_numbered.append((header_with_rep, genome_length, all_genes))

print(f"Processed {len(genomes_data_numbered)} genomes")

binary_out = f"{args.out_prefix}_presence_binary_data_fromTSV_{formatted_evalue}"
with open(binary_out, "w") as f:
    for item in genomes_data_numbered:
        f.write(str(item) + "\n")

# ----- Reflect synteny around cog_of_interest -----
print("Creating reflect synteny...")
import ast

# Definition of cog_of_interest-------------------------------------------
print(f"cog_of_interest: {cog_of_interest}")
#locustag_of_interest = locus_tag_cog_dic.get(cog_of_interest, None)
#print(locustag_of_interest)


def reflect_synteny(genes):
    goi = next((g for g in genes if any(locus_tag.lower() in g[4].lower() for locus_tag in locus_tags_cog_list)), None)
    if not goi:
        return genes
    print(f"Reflect synteny{goi}")
    goi_start, goi_end, goi_strand = goi[0], goi[1], goi[2]
    goi_center = (goi_start + goi_end) // 2
    if goi_strand == 1:
        return genes
    reflected_genes = []
    for g in genes:
        start, end, strand, name, tag = g
        new_start = 2 * goi_center - end
        new_end = 2 * goi_center - start
        new_strand = -strand
        reflected_genes.append((min(new_start, new_end), max(new_start, new_end), new_strand, name, tag))
    reflected_genes.sort(key=lambda x: x[0])
    return reflected_genes

input_path = binary_out
oriented_out = f"{args.out_prefix}_genomic_context_data_{cog_of_interest}_{formatted_evalue}_oriented"
with open(input_path, "r") as infile, open(oriented_out, "w") as outfile:
    for line in infile:
        data = ast.literal_eval(line.strip())
        header, genome_length, genes = data
        genes_oriented = reflect_synteny(genes)
        outfile.write(f"{(header, genome_length, genes_oriented)}\n")

# ----- Extract context window around cog_of_interest -----
print("Creating context gene...")

def context_goi(genes):
    #goi = next((g for g in genes if cog_of_interest in g[3].lower()), None)
    #goi = next((g for g in genes if locustag_of_interest in g[3].lower()), None)
    goi = next((g for g in genes if any(locus_tag.lower() in g[4].lower() for locus_tag in locus_tags_cog_list)), None)
   
    print(f"context goi {goi}")
    if not goi:
        return genes
    goi_start, goi_end, _ = goi[0], goi[1], goi[2]
    goi_center = (goi_start + goi_end) // 2
    goi_context_upstream = goi_center - 60000
    goi_context_downstream = goi_center + 60000
    context_genes = []
    for g in genes:
        start, end, strand, name, tag = g
        if goi_context_upstream < start < goi_context_downstream:
            context_genes.append(g)
    context_genes.sort(key=lambda x: x[0])
    return context_genes

sector_out = f"{args.out_prefix}_genomic_context_data_{cog_of_interest}_{formatted_evalue}_sector"
with open(oriented_out, "r") as infile, open(sector_out, "w") as outfile:
    for line in infile:
        data = ast.literal_eval(line.strip())
        header, genome_length, genes = data
        genes_oriented = context_goi(genes)
        outfile.write(f"{(header, genome_length, genes_oriented)}\n")