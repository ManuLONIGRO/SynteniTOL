#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates an iTOL DATASET_ARROWS synteny file centered on the gene linked to
each protein_id present in the input FASTA.

Unlike createFiles2Synteny.py (which uses density/centrality scoring to pick
the best GOI per organism), this script maps each FASTA protein_id directly
to the CDS protein_id qualifiers in the .gbff files, extracts the genomic
context around that gene, reflects the synteny to a common orientation, and
writes the iTOL dataset in a single self-contained step. Protein_ids for
which no exact matching CDS is found are skipped with a warning.
"""

import os
import argparse
from typing import Dict, List, Tuple

from color_palettes import categorical_colors

import pandas as pd
from Bio import SeqIO


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Generate iTOL synteny ARROWS centered on FASTA protein_ids."
)
parser.add_argument(
    "--input_fasta", type=str, required=True,
    help="Path to input FASTA file (protein_ids parsed from headers).",
)
parser.add_argument(
    "--input_dir", type=str, required=True,
    help="Directory containing .gbff files.",
)
parser.add_argument(
    "--dataframe", type=str, required=True,
    help="Path to dataframe.tsv (for locus-to-gene mapping and COG coloring).",
)
parser.add_argument(
    "--evalue", type=float, default=1e-10,
    help="E-value threshold for filtering dataframe hits (default: 1e-10).",
)
parser.add_argument(
    "--cogs", type=str, default=None,
    help="Comma-separated COG/gene list for coloring, e.g. COG1152,COG1795.",
)
parser.add_argument(
    "--gene_list", type=str, nargs="*", default=None,
    help="Space-separated gene/COG list for coloring (alternative to --cogs).",
)
parser.add_argument(
    "--color_by_group", type=str, default=None,
    help="Comma-separated groups; hyphen joins members sharing a color, "
    "e.g. COG1152-COG1229,COG1795.",
)
parser.add_argument(
    "--itol_synteny_file", type=str, required=True,
    help="Output path for iTOL DATASET_ARROWS file.",
)
args = parser.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def parse_cogs(cogs_str: str) -> List[str]:
    if not cogs_str:
        return []
    return [c.strip() for c in cogs_str.split(",") if c.strip()]


def parse_color_groups(spec: str) -> Dict[str, List[str]]:
    groups: Dict[str, List[str]] = {}
    if not spec:
        return groups
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        members = [p.strip() for p in token.split("-") if p.strip()]
        if members:
            groups[token] = members
    return groups


def assign_colors(
    genes: List[str],
    color_groups: Dict[str, List[str]],
) -> Tuple[Dict[str, str], List[str]]:
    keys: List[str] = []
    member_to_group: Dict[str, str] = {}
    for group_key, members in color_groups.items():
        keys.append(group_key)
        for m in members:
            member_to_group[m] = group_key
    for g in genes:
        if g not in member_to_group and g not in keys:
            keys.append(g)
    palette = categorical_colors(len(keys))
    key_to_color = {k: palette[i] for i, k in enumerate(keys)}
    gene_to_color = {
        g: key_to_color.get(member_to_group.get(g, g), "#000000")
        for g in genes
    }
    ordered_colors = [gene_to_color[g] for g in genes]
    return gene_to_color, ordered_colors


def clean_metadata(value: str) -> str:
    return (
        value.replace(" ", "_")
        .replace(";", "_")
        .replace(":", "_")
        .replace("=", "_")
        .replace("(", "_")
        .replace(")", "_")
    )


# ---------------------------------------------------------------------------
# STEP 1: Parse protein_ids from FASTA
# ---------------------------------------------------------------------------
print(f"Reading protein_ids from {args.input_fasta}")
fasta_protein_ids: List[str] = []
for record in SeqIO.parse(args.input_fasta, "fasta"):
    header_token = record.id.split(":")[0].split()[0].strip()
    if header_token:
        fasta_protein_ids.append(header_token)

print(f"  Found {len(fasta_protein_ids)} protein_id(s) in FASTA")
if not fasta_protein_ids:
    raise SystemExit("No protein_ids found in the input FASTA file.")


# ---------------------------------------------------------------------------
# STEP 2: Read dataframe, build locus-to-gene lookup
# ---------------------------------------------------------------------------
print(f"Reading dataframe from {args.dataframe}")
df = pd.read_csv(args.dataframe, sep="\t", dtype={"protein_id": str})
df_filtered = df[df["acc"] >= 0.60].copy()

locus_to_gene: Dict[str, str] = dict(
    zip(df_filtered["locus_tag"], df_filtered["gene"])
)

# Build gene list for coloring
gene_list: List[str] = []
if args.gene_list:
    gene_list = [name.strip() for name in args.gene_list if str(name).strip()]
if not gene_list and args.cogs:
    gene_list = parse_cogs(args.cogs)

color_groups: Dict[str, List[str]] = {}
if args.color_by_group:
    color_groups = parse_color_groups(args.color_by_group)
genes_color_dic, _gen_color = assign_colors(gene_list, color_groups)
print(f"  Gene list for coloring: {gene_list}")
print(f"  Color groups: {color_groups}")


# ---------------------------------------------------------------------------
# STEP 3: Match protein_ids to CDS protein_id qualifiers in the .gbff files
# ---------------------------------------------------------------------------
print(f"Matching protein_ids to CDS protein_id qualifiers in {args.input_dir}")
gbff_files = [f for f in os.listdir(args.input_dir) if f.endswith(".gbff")]
fasta_protein_set = set(fasta_protein_ids)

# protein_id → (contig_genes, center_locus_tag, taxid, species, strain, isolate)
matched: Dict[str, Tuple[List[Tuple], str, str, str, str, str]] = {}

for gbff_name in gbff_files:
    gbff_path = os.path.join(args.input_dir, gbff_name)
    records = list(SeqIO.parse(gbff_path, "genbank"))
    if not records:
        print(f"  WARNING: could not parse {gbff_name}, skipped")
        continue

    for record in records:
        contig_genes: List[Tuple[int, int, int, str, str]] = []
        taxid = "NA"
        species = "NA"
        strain = "NA"
        isolate = "NA"
        pid_to_locus: Dict[str, str] = {}

        for feature in record.features:
            if feature.type == "source":
                db_xrefs = feature.qualifiers.get("db_xref", [])
                isolate = clean_metadata(
                    feature.qualifiers.get("isolate", ["NA"])[0]
                )
                strain = clean_metadata(
                    feature.qualifiers.get("strain", ["NA"])[0]
                )
                species = clean_metadata(
                    feature.qualifiers.get("organism", ["NA"])[0]
                )
                for xref in db_xrefs:
                    if xref.startswith("taxon:"):
                        taxid = xref.split(":")[1]
                        break

            if feature.type == "CDS":
                lt = feature.qualifiers.get("locus_tag", ["NA"])[0]
                p_id = feature.qualifiers.get("protein_id", ["NA"])[0].strip()
                if not lt:
                    continue
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = int(feature.location.strand)
                gene_name = locus_to_gene.get(lt, lt)
                contig_genes.append((start, end, strand, gene_name, lt))

                if p_id in fasta_protein_set:
                    pid_to_locus[p_id] = lt

        for pid, center_locus_tag in pid_to_locus.items():
            if (
                pid not in matched
                or len(contig_genes) > len(matched[pid][0])
            ):
                matched[pid] = (
                    list(contig_genes),
                    center_locus_tag,
                    taxid,
                    species,
                    strain,
                    isolate,
                )

not_found = [pid for pid in fasta_protein_ids if pid not in matched]
if not_found:
    print(
        f"  WARNING: {len(not_found)} protein_id(s) not found as CDS in "
        f"any .gbff and will be skipped: {', '.join(not_found)}"
    )
print(f"  Matched {len(matched)} / {len(fasta_protein_ids)} protein_ids")


# ---------------------------------------------------------------------------
# STEP 4: For each matched protein_id, center the context on its CDS
# ---------------------------------------------------------------------------
results: List[Tuple[str, int, List[Tuple], str]] = []

for pid in fasta_protein_ids:
    hit = matched.get(pid)
    if hit is None:
        print(f"  WARNING: protein_id {pid} not found as CDS, skipped")
        continue

    contig_genes, center_locus_tag, taxid, species, strain, isolate = hit
    if not contig_genes:
        print(f"  WARNING: no contig genes found for {pid}")
        continue

    # -- reflect synteny --
    goi = next(
        (g for g in contig_genes if g[4] == center_locus_tag), None
    )
    if goi is None:
        print(f"  WARNING: center locus_tag {center_locus_tag} not in contig genes")
        continue

    genes = list(contig_genes)
    goi_start, goi_end, goi_strand = goi[0], goi[1], goi[2]
    if goi_strand == -1:
        goi_center = (goi_start + goi_end) // 2
        reflected = []
        for g in genes:
            s, e, st, nm, tg = g
            ns = 2 * goi_center - e
            ne = 2 * goi_center - s
            reflected.append((min(ns, ne), max(ns, ne), -st, nm, tg))
        genes = sorted(reflected, key=lambda x: x[0])

    # -- context window (50 kb each side of the gene = 100 kb total) --
    goi_ref = next((g for g in genes if g[4] == center_locus_tag), None)
    if goi_ref is None:
        continue
    ctx_start = goi_ref[0] - 50000
    ctx_end = goi_ref[1] + 50000
    context_genes = [g for g in genes if ctx_start <= g[0] <= ctx_end]
    context_genes.sort(key=lambda x: x[0])

    if not context_genes:
        print(f"  WARNING: empty context for {pid}")
        continue

    header = f"{taxid}|{species}|{strain}|{isolate}"
    results.append((header, 120000, context_genes, center_locus_tag))
    print(f"  OK: {pid} -> {header} ({len(context_genes)} genes in context)")


# ---------------------------------------------------------------------------
# STEP 5: Generate iTOL DATASET_ARROWS file
# ---------------------------------------------------------------------------
label = f"Synteny_{fasta_protein_ids[0]}"
out_path = args.itol_synteny_file

if not results:
    print(
        "\nWARNING: 0 protein_id(s) could be located in the input genomes. "
        "An empty iTOL synteny file (header only) will be written so the "
        "pipeline does not fail."
    )

with open(out_path, "w") as fout:
    fout.write(
        f"DATASET_ARROWS\n"
        f"SEPARATOR COMMA\n"
        f"DATASET_LABEL,{label}\n"
        f"COLOR,#0000ff\n"
        f"WIDTH,3000\n"
        f"BORDER_WIDTH,0.5\n"
        f"\n"
        f"DATA\n"
    )

    if results:
        seen_headers: Dict[str, int] = {}
        for header, window, context_genes, center_lt in results:
            seen_headers[header] = seen_headers.get(header, 0) + 1
            org_line = f"{header}|{seen_headers[header]}"

            center_goi = next(
                (g for g in context_genes if g[4] == center_lt), None
            )
            if center_goi is None:
                continue

            start_goi = center_goi[0]
            entries: List[str] = []
            for g in context_genes:
                start, end, strand, gene_name, locus_tag = g
                color = genes_color_dic.get(gene_name, "#FFFFFF")
                start_ref = start - start_goi + 50000
                end_ref = end - start_goi + 50000
                gene_size = abs(end_ref - start_ref)
                head_arrow = (
                    300 if gene_size > 310
                    else 200 if gene_size > 210
                    else 100 if gene_size > 110
                    else 30
                )
                if strand == -1:
                    hl, hr = str(head_arrow), "0"
                else:
                    hl, hr = "0", str(head_arrow)
                entries.append(
                    f"{start_ref}|{end_ref}|{color}|{color}|"
                    f"{gene_name}|#000000|1|0.5|{hl}|{hr}"
                )

            fout.write(f"{org_line},{window},{','.join(entries)}\n")

if results:
    print(f"\nWrote {len(results)} organism(s) to {out_path}")
else:
    print(f"\nWrote empty iTOL synteny file (0 organisms) to {out_path}")
