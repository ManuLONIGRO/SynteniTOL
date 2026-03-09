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

parser = argparse.ArgumentParser()
parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold for filtering")
parser.add_argument("--input_dir", type=str, required=True, help="Input directory containing .gbff files")
parser.add_argument("--dataframe", type=str, required=True, help="Path to the dataframe TSV file")
parser.add_argument("--out_prefix", type=str, required=True, help="Output prefix for generated files")
parser.add_argument("--gene_list", type=str, nargs='*', default=None, help="List of genes/COGs of interest (space-separated)")
parser.add_argument("--cogs", type=str, default=None, help="Comma-separated COG list, e.g., COG1152,COG1795")
parser.add_argument("--protein_to_organism_map_tsv", type=str, required=True, help="Path to protein_to_organism_map.tsv file")
parser.add_argument("--out_tsv", type=str, required=True, help="Path to output TSV file with species info")
parser.add_argument("--out_candidates_tsv", type=str, required=False, help="Path to output candidates TSV file")
parser.add_argument("--out_best_goi_tsv", type=str, required=False, help="Path to output best GOI per organism TSV file")
args = parser.parse_args()

choose_evalue = args.evalue
formatted_evalue = f"{choose_evalue:.0e}"  # scientific notation for filenames

# Load dataframe and build lookups
print(f"Reading dataframe from {args.dataframe}")

df = pd.read_csv(args.dataframe, sep="\t", dtype={'protein_id': str})
# Filter by coverage
df_filtered = df[df["coverage"] > 0.60]

# Build lookups
locus_to_gene: Dict[str, str] = dict(zip(df_filtered["locus_tag"], df_filtered["gene"]))
#print(locus_to_gene)
#dic_taxid: Dict[str, str] = dict(zip(df_filtered["locus_tag"], df_filtered["taxid"]))
#dic_species: Dict[str, str] = dict(zip(df_filtered["locus_tag"], df_filtered["specie"]))
#dic_isolate: Dict[str, str] = dict(zip(df_filtered["locus_tag"], df_filtered["isolate"]))

print(f"Creating presence_binary_data_fromTSV_{formatted_evalue}...")

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
#---------------------------------------------------------------------------------------------------------------------------
# Select the best hit per organism for the cog_of_interest NEW LINES
df_cog= df_filtered[
    (df_filtered["gene"] ==  cog_of_interest) &
    (df_filtered["evalue"] <= choose_evalue)].copy()
df_cog["organism"] = df_cog["locus_tag"].str.split("_").str[0]
# df_cog["position_in_genome"] = df_cog["locus_tag"].str.split("_").str[1].astype(int)
#%%
# # Function to pick best locus_tag by density of genes in the context window of 600
# def best_hit_per_organism_by_prefix(
#     df,
#     cog_of_interest=cog_of_interest,
#     window=300,
#     locus_col="locus_tag",
#     cog_col="gene",
# ):
#     df = df.copy()

#     # 1) organismo = prefijo antes del primer "_"
#     df["organism"] = df[locus_col].astype(str).str.split("_", n=1).str[0]

#     # 2) número del locus_tag (últimos dígitos)
#     # df["locus_num"] = df[locus_col].astype(str).str.extract(r"(\d+)$")[0]
#     # df = df.dropna(subset=["locus_num"])
#     # df["locus_num"] = df["locus_num"].astype(int)
#     df["locus_num"] = (
#     df[locus_col]
#     .astype(str)
#     .str.extract(r"(\d+)(?!.*\d)", expand=False)  # último grupo numérico del string
# )
#     df = df.dropna(subset=["locus_num"])
#     df["locus_num"] = df["locus_num"].astype(int)

#     # 3) filtrar candidatos (solo filas del COG de interés)
#     candidates = df[df[cog_col] == cog_of_interest].copy()

#     if candidates.empty:
#         return pd.DataFrame()

#     # 4) contar cuántas proteínas hay en la ventana ±window para cada candidato
#     counts = []
#     for idx, row in candidates.iterrows():
#         org = row["organism"]
#         center = row["locus_num"]
#         lo, hi = center - window, center + window

#         sub = df[df["organism"] == org]
#         n_total = ((sub["locus_num"] >= lo) & (sub["locus_num"] <= hi)).sum()

#         #counts.append((idx, n_total))
#         counts.append((n_total))

#     #candidates["window_count_total"] = pd.Series(dict(counts))
#     candidates["window_count_total"] = counts

#     # 5) elegir el mejor por organismo (si hay empates, queda el primero)
#     best = (
#         candidates.sort_values(["organism", "window_count_total"], ascending=[True, False])
#         .groupby("organism", as_index=False)
#         .head(1)
#         .reset_index(drop=True)
#     )

#     return best, candidates.sort_values(["organism", "window_count_total"], ascending=[True, False])

# # Use de function to pick the best locust_tag for the cog of interest
# best_goi_dataframe, candidates_dataframe = best_hit_per_organism_by_prefix(df_filtered, cog_of_interest, 300)
#%%
def best_hit_per_organism_by_prefix(
    df,
    cog_of_interest=cog_of_interest,
    gene_list=None,
    window=300,
    locus_col="locus_tag",
    cog_col="gene",
    evalue_threshold=1e-60,
    coverage_threshold=0.90
):
    df = df.copy()
    if gene_list is None:
        gene_list = [cog_of_interest]

    # 1) Define organism and extract numeric locus part
    df["organism"] = df[locus_col].astype(str).str.split("_", n=1).str[0]
    df["locus_num"] = (
        df[locus_col]
        .astype(str)
        .str.extract(r"(\d+)(?!.*\d)", expand=False) 
    )
    df = df.dropna(subset=["locus_num"])
    df["locus_num"] = df["locus_num"].astype(int)

    # 2) Filter candidates (rows matching the specific COG of interest)
    candidates = df[df[cog_col] == cog_of_interest].copy()
    if candidates.empty:
        return pd.DataFrame(), pd.DataFrame()

    # 3) Calculate Density and Centrality
    # Identify High Confidence hits (Tiering)
    # True (1) is sorted higher than False (0) when using ascending=False
    candidates["is_high_confidence"] = candidates["evalue"] <= evalue_threshold
    density_counts = []
    centrality_scores = []

    for idx, row in candidates.iterrows():
        org = row["organism"]
        center = row["locus_num"]
        lo, hi = center - window, center + window

        # Get all "target" genes in this organism's window
        sub = df[(df["organism"] == org) & (df["gene"].isin(gene_list))]
        neighbors = sub[(sub["locus_num"] >= lo) & (sub["locus_num"] <= hi)]
        
        # Density: count of target genes
        n_total = len(neighbors)
        density_counts.append(n_total)

        # Centrality: Mean absolute distance to neighbors 
        # (Lower is better/more centered)
        if n_total > 1:
            mean_dist = (neighbors["locus_num"] - center).abs().mean()
        else:
            mean_dist = 0 # Solo gene is technically centered
        centrality_scores.append(mean_dist)

    candidates["window_count_total"] = density_counts
    candidates["centrality_score"] = centrality_scores

    # 4) The "Context-First" Sorting Strategy
    # - is_high_confidence: Descending (High confidence hits first)
    # - window_count_total: Descending (Most neighbors next)
    # - centrality_score: Ascending (Most centered)
    # - evalue: Ascending (Final tie-break match quality)
    # - locus_num: Ascending (Deterministic fallback)
    
    # 5) identify good coverage hits
    candidates["is_good_coverage"] = candidates["coverage"] >= coverage_threshold

    sorted_candidates = candidates.sort_values(
        by=["organism", "is_high_confidence", "is_good_coverage", "window_count_total", "centrality_score", "evalue", "locus_num"],
        ascending=[True, False, False ,False, True, True, True]
    )

    # 6) Choose the best one per organism
    best = sorted_candidates.groupby("organism", as_index=False).head(1).reset_index(drop=True)

    return best, sorted_candidates

# Usage remains the same:
best_goi_dataframe, candidates_dataframe = best_hit_per_organism_by_prefix(
    df_filtered, 
    cog_of_interest, 
    gene_list=gene_list, # Passes the whole list for density/centrality
    window=300
)
#%%
#output candidates in a new tsv
candidates_flag = args.out_candidates_tsv
if candidates_flag:
    candidates_dataframe.to_csv(candidates_flag, sep="\t", index=False)

best_goi_dataframe_tsv_flag = args.out_best_goi_tsv
if best_goi_dataframe_tsv_flag:
    best_goi_dataframe.to_csv(best_goi_dataframe_tsv_flag, sep="\t", index=False)


#----------------------------------------------------------------------------------------------------------------------------
# Define organism from locus_tag prefix 
df_cog["organism"] = df_cog["locus_tag"].str.split("_").str[0]
# Sort by evalue ascending
df_cog = df_cog.sort_values(by="evalue", ascending=True)

# One hit per organism, keep the  first
df_best_per_org = df_cog.drop_duplicates(subset=["organism"], keep="first")

# These are the protin ids and locus_tags of the best result
#protein_interest_list = df_best_per_org["protein_id"].tolist() 03/0
protein_interest_list = best_goi_dataframe["protein_id"].tolist()
protein_interest_set = set(protein_interest_list)
#locus_tag_goi = df_best_per_org["locus_tag"].tolist()
locus_tag_goi = best_goi_dataframe["locus_tag"].tolist()

#------------------------------------------------------------------------------------
# #Ordenar por evalue ascendente
# df_cog = df_cog.sort_values(by="evalue", ascending=True)
# # quedarme con 1 hit por organismo
# df_best_per_org = df_cog.drop_duplicates(subset=["organism"], keep="first")
# protein_interest_list = df_best_per_org["protein_id"].tolist()
# protein_interest_set = set(protein_interest_list)
# locus_tag_goi = df_best_per_org["locus_tag"].tolist() #new line

#------------------------------------------------------------------------------------
#protein_interest_list = df_filtered[df_filtered["gene"] == cog_of_interest]["protein_id"].unique().tolist()
#protein_interest_set = set([pid for pid in protein_interest_list])
print(f"Proteins of interest for {cog_of_interest}: {protein_interest_set}")
#this works
#locus_tags_cog_list = df_filtered[["gene"] == cog_of_interest]["locus_tag"].unique().tolist() #works
locus_tags_cog_list = df_filtered[df_filtered["gene"] == cog_of_interest]["locus_tag"].unique().tolist()
locus_goi_set = locus_tags_cog_list #works

#print(f"locus_tag cog list {locus_goi_set}")
#locus_tags_cog_list = df_filtered[df_filtered["gene"] == cog_of_interest]["locus_tag"].unique().tolist()
# locus_goi_set = set([lt for lt in locus_tags_cog_list])

# Rows for the new tsv with species info
rows_to_tsv = []
# Load protein_to_organism_map.tsv to build dict protein: organism,isolate

def process_gbff(input_file: str, locus_to_gene: Dict[str, str], locus_goi_set: set):

    records = list(SeqIO.parse(input_file, "genbank"))
    if not records:
        return None

    best_contig_genes: List[Tuple[str, str, str, str, str]] = []
    best_contig_len = 0
    taxid = None
    # isolate = "NA"

    # choose only the contig that contains the GOI; if multiple, pick the one with more genes
    for record in records:
        contig_genes: List[Tuple[str, str, str, str, str]] = []
        contains_goi = False
        local_taxid = "NA"
        local_isolate = "NA"
        local_strain = "NA"
        local_species = "NA"
        local_protein_id = "NA"

        for feature in record.features:
            if feature.type == "source":
                db_xrefs = feature.qualifiers.get("db_xref", [])
                local_isolate = feature.qualifiers.get("isolate", ["NA"])[0].replace(" ", "_").replace(";","").replace(":","_").replace("(","_").replace(")","_")
                local_strain = feature.qualifiers.get("strain", ["NA"])[0].replace(" ", "_").replace(";","_").replace(":","_").replace("(","_").replace(")","_")
                local_species = feature.qualifiers.get("organism", ["NA"])[0].replace(" ", "_").replace(";","").replace(":","_").replace("(","_").replace(")","_")
                for xref in db_xrefs:
                    if xref.startswith("taxon:"):
                        local_taxid = xref.split(":")[1]
                        break
            if feature.type == "CDS":
                locus_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]
                # local_protein_id = feature.qualifiers.get("protein_id", ["NA"])[0].split(".")[0].strip()
                local_protein_id = feature.qualifiers.get("protein_id", ["NA"])[0].strip()
                if not locus_tag:
                    continue
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = int(feature.location.strand)
                gene_name = locus_to_gene.get(locus_tag, locus_tag)
                contig_genes.append((start, end, strand, gene_name, locus_tag))
                #if locus_tag in locus_goi_set:
                if locus_tag in locus_tag_goi: #new line
                    contains_goi = True
                if local_protein_id in protein_interest_set:
                    print(f"Found local protein id {local_protein_id} in protein_interest_set")
                    rows_to_tsv.append({
                        "protein_id": local_protein_id,
                        "local_taxid": local_taxid,
                        "locus_tag": locus_tag,
                        "local_species": local_species,
                        "local_strain": local_strain,
                        "local_isolate": local_isolate
                    })
                    # to_map = f"{local_protein_id}\t{local_taxid}\t{locus_tag}\t{local_species}\t{local_strain}\t{local_isolate}"
    
        if contains_goi and len(contig_genes) > len(best_contig_genes):
            best_contig_genes = contig_genes
            best_contig_len = len(record.seq)
            taxid = local_taxid
            
            # isolate = local_isolate
        # if local_protein_id:
        #     header = local_protein_id
        # else:
        #     header = "Unknown|Unknown|NA"
    

    if not best_contig_genes or not taxid:
        return None
    if local_species:
        header = f"{local_taxid}|{local_species}|{local_strain}|{local_isolate}"
    else:
        header = "Unknown|Unknown|Unknown|Unknown|NA"
    # header = f"{species_taxid[0]}|{species_name}|{isolate}"
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


# def reflect_synteny(genes):
#     goi = next((g for g in genes if any(locus_tag in g[4] for locus_tag in locus_tags_cog_list)), None)
#     if not goi:
#         return genes
#     print(f"Reflect synteny{goi}")
#     goi_start, goi_end, goi_strand = goi[0], goi[1], goi[2]
#     goi_center = (goi_start + goi_end) // 2
#     if goi_strand == 1:
#         return genes
#     reflected_genes = []
#     for g in genes:
#         start, end, strand, name, tag = g
#         new_start = 2 * goi_center - end
#         new_end = 2 * goi_center - start
#         new_strand = -strand
#         reflected_genes.append((min(new_start, new_end), max(new_start, new_end), new_strand, name, tag))
#     reflected_genes.sort(key=lambda x: x[0])
#     return reflected_genes

def reflect_synteny(genes):
    # use only the selected best-per-organism GOI locus_tags
    goi = next((g for g in genes if g[4] in locus_tag_goi), None)
    if not goi:
        return genes
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

# def context_goi(genes):
#     goi = next((g for g in genes if any(locus_tag in g[4] for locus_tag in locus_tags_cog_list)), None)
    
#     print(f"context goi {goi}")
    
#     if not goi:
#         return genes
#     goi_start, goi_end, _ = goi[0], goi[1], goi[2]
#     goi_center = (goi_start + goi_end) // 2
#     goi_context_upstream = goi_center - 60000
#     goi_context_downstream = goi_center + 60000
#     context_genes = []
#     for g in genes:
#         start, end, strand, name, tag = g
#         if goi_context_upstream < start < goi_context_downstream:
#             context_genes.append(g)
#     context_genes.sort(key=lambda x: x[0])
#     return context_genes

def context_goi(genes):
    # again, only the selected GOI locus_tags
    goi = next((g for g in genes if g[4] in locus_tag_goi), None)

    print(f"context goi {goi}")

    if not goi:
        return genes
    goi_start, goi_end, _ = goi[0], goi[1], goi[2]
    goi_center = (goi_start + goi_end) // 2
    goi_context_upstream = goi_center - 100000
    goi_context_downstream = goi_center + 100000
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

# Load protein_to_organism_map.tsv
map_tsv = args.protein_to_organism_map_tsv
df_organism_map = pd.read_csv(map_tsv, sep="\t").rename(columns={'protein': 'protein_id'})

# Create the new tsv with species info
df_new = pd.DataFrame(rows_to_tsv)
df_final = df_organism_map.merge(df_new, on="protein_id", how="left")
out_tsv = args.out_tsv
df_final.to_csv(out_tsv, sep="\t", index=False)