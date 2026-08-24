#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 23 10:51:44 2025

Updated for NEXTFLOW

- Reads synteny inputs produced by createFiles2Synteny.py and organism gene presence from itol_binary.py
- Dynamic gene list from --cogs (comma-separated) or --gene_list
- Optional --color_by_group: e.g. COG1152-COG1229,COG1795 groups colors by token
- Generates contrasting colors automatically

Header format expected in data files: taxid|species_name|isolate|n_repetitions
"""

import ast
import argparse
from typing import Dict, List, Tuple
import pandas as pd

from color_palettes import categorical_colors


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold for filtering")
    # Either pass gene names directly or COG list; at least one must be provided
    parser.add_argument("--gene_list", type=str, nargs='*', default=None, help="List of genes/COGs of interest (space-separated)")
    parser.add_argument("--cogs", type=str, default=None, help="Comma-separated COG list, e.g., COG1152,COG1795")
    parser.add_argument("--reference_gene", type=str, default=None, help="Gene used as synteny reference in labels/logs (defaults to first gene of the list)")
    parser.add_argument("--color_by_group", type=str, default=None, help="Comma-separated groups; hyphen joins members, e.g. COG1152-COG1229,COG1795")
    parser.add_argument("--genomic_context_data", type=str, required=True, help="Path to genomic context data file")
    parser.add_argument("--itol_synteny_file", type=str, required=True, help="Output path for itol synteny oriented file")
    parser.add_argument("--synteny_log", type=str, required=True, help="Output path for synteny log file")
    parser.add_argument("--organism_genes_tsv", type=str, required=True, help="Path to organism gene presence TSV from itol_binary.py")
    parser.add_argument("--itol_profiling_file", type=str, required=True, help="Output path for itol binary data file")
    parser.add_argument("--best_goi_tsv", type=str, required=True, help="Path to best GOI per organism TSV file")
    return parser.parse_args()


def parse_cogs(cogs_str: str) -> List[str]:
    if not cogs_str:
        return []
    return [c.strip() for c in cogs_str.split(',') if c.strip()]

def parse_color_groups(spec: str) -> Dict[str, List[str]]:
    groups: Dict[str, List[str]] = {}
    if not spec:
        return groups
    for token in spec.split(','):
        token = token.strip()
        if not token:
            continue
        members = [p.strip() for p in token.split('-') if p.strip()]
        if members:
            groups[token] = members
    return groups


def assign_colors(genes: List[str], color_groups: Dict[str, List[str]]) -> Tuple[Dict[str, str], List[str]]:
    # Create keys for coloring: group tokens or individual genes
    keys: List[str] = []
    member_to_group: Dict[str, str] = {}
    for group_key, members in color_groups.items():
        keys.append(group_key)
        for m in members:
            member_to_group[m] = group_key
    # Add standalone genes not in any group
    for g in genes:
        if g not in member_to_group and g not in keys:
            keys.append(g)
    # Generate palette
    palette = categorical_colors(len(keys))
    key_to_color = {k: palette[i] for i, k in enumerate(keys)}
    # Map each gene to a color via its group (if present) else itself
    gene_to_color = {g: key_to_color.get(member_to_group.get(g, g), "#000000") for g in genes}
    # The FIELD_COLORS order should match genes order
    ordered_colors = [gene_to_color[g] for g in genes]
    return gene_to_color, ordered_colors


# -------------- Main --------------
args = parse_args()

evalue_limit = args.evalue
formatted_evalue = f"{evalue_limit:.0e}"

# Build gene list
gene_list: List[str] = []
if args.gene_list:
    gene_list = list(args.gene_list)
if not gene_list and args.cogs:
    gene_list = parse_cogs(args.cogs)
if not gene_list:
    raise SystemExit("You must provide --gene_list or --cogs")

gene_of_interest = args.reference_gene if args.reference_gene else gene_list[0]

# Parse color grouping
color_groups = parse_color_groups(args.color_by_group) if args.color_by_group else {}

# Assign colors
genes_color_dic, gen_color = assign_colors(gene_list, color_groups)

# Inputs/outputs
genomes_itol = args.genomic_context_data
archivo_sintenia = args.itol_synteny_file
archivo_log = args.synteny_log
archivo_binario = args.organism_genes_tsv
itol_binario = args.itol_profiling_file
headers_file = f"organismos_binario_{gene_of_interest}_{formatted_evalue}"

# ---------- Load context data ----------
with open(genomes_itol, "r") as f:
    lineas = f.readlines()

datos = []
for linea in lineas:
    parsed_line = ast.literal_eval(linea.strip())
    organism_name = parsed_line[0]
    genome_length = parsed_line[1]
    genes = parsed_line[2]
    datos.append((organism_name, genome_length, genes))

# ----------- Load best_goi.tsv ----------------
best_goi_tsv = args.best_goi_tsv
best_goi_df = pd.read_csv(best_goi_tsv, sep="\t")
best_goi_locustags_list = best_goi_df["locus_tag"].tolist()

# ---------- Synteny (iTOL domains) ----------
with open(archivo_sintenia, "w") as f, open(archivo_log, "w") as log_file:
    f.write(f"""DATASET_ARROWS
SEPARATOR COMMA
DATASET_LABEL,Synteny_{gene_of_interest}
COLOR,#0000ff
WIDTH,3000
BORDER_WIDTH,0.5

DATA
""")
    for item in datos:
        organism_name = item[0]
        genome_length = item[1]
        genes = item[2]

        gene_context = []
        context_inicio = None
        context_fin = None
        start_goi = None

        # Find GOI to set context window
        for gene in genes:
            start, end, strand, gene_name, locus_tag = gene
            if locus_tag in best_goi_locustags_list:
                print(f"Found {locus_tag} in {organism_name} as GOI at {start}-{end} (strand {strand})")
                start_goi = start
                context_inicio = start - 50000
                context_fin = end + 50000
                break

        if start_goi is None:
            # No GOI in this organism; skip
            continue

        # Collect genes mapped to window with coloring
        for gene in genes:
            start, end, strand, gene_name, locus_tag = gene
            color = genes_color_dic.get(gene_name, "#FFFFFF")
            if context_inicio <= start <= context_fin:
                start_reference = start - start_goi + 50000
                end_reference = end - start_goi + 50000
                head_arrow = 300 if abs(end_reference - start_reference) > 310 else 200 if abs(end_reference - start_reference) > 210 else 100 if abs(end_reference - start_reference) > 110 else 30  # Short genes get smaller heads avoiding an error in iTOL
                # iTOL arrow direction is controlled by which side has a non-zero head width:
                # - head on the LEFT  => arrow points LEFT  (strand -1)
                # - head on the RIGHT => arrow points RIGHT (strand +1)
                direction, head_width_left, head_width_right = ("left", head_arrow, "0") if strand == -1 else ("right", "0", head_arrow)
                gene_entry_context = f"{start_reference}|{end_reference}|{color}|{color}|{gene_name}|#000000|1|0.5|{head_width_left}|{head_width_right}"
                gene_context.append(gene_entry_context)

        log_file.write(f"Context of {gene_of_interest} of {organism_name}\n{gene_context}\n\n")
        f.write(f"{organism_name},120000,{','.join(gene_context)}\n")

# ---------- Binary profiling ----------
organism_gene_presence: Dict[str, Dict[str, str]] = {}
organism_names: List[str] = []

presence_df = pd.read_csv(archivo_binario, sep="\t", dtype=str)
for _, row in presence_df.iterrows():
    organism_name = row["organism"]
    genes_present = {g.strip() for g in str(row["genes"]).split(",") if g.strip()}
    organism_names.append(organism_name)
    organism_gene_presence[organism_name] = {
        gene: ("1" if gene in genes_present else "-1") for gene in gene_list
    }

n_fields = ",1" * len(gene_list)

with open(itol_binario, "w") as f3:
    f3.write(
        f"DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,Gene\nCOLOR,#ff0000\nFIELD_SHAPES{n_fields}\nFIELD_LABELS,{','.join(gene_list)}\nFIELD_COLORS,{','.join(gen_color)}\nDATA\n"
    )
    for organism, gene_presence in organism_gene_presence.items():
        presence_list = [str(gene_presence[gene]) for gene in gene_list]
        f3.write(f"{organism},{','.join(presence_list)}\n")

with open(headers_file, "w") as f4:
    for name in organism_names:
        f4.write(f"{name}\n")