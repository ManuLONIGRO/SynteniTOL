#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates an iTOL DATASET_ARROWS synteny file centered on the gene linked to
each protein_id present in the input FASTA.

Unlike createFiles2Synteny.py (which uses density/centrality scoring to pick
the best GOI per organism), this script maps each FASTA protein_id directly
to its .gbff file via the filename, extracts the genomic context around that
gene, reflects the synteny to a common orientation, and writes the iTOL
dataset in a single self-contained step.
"""

import os
import argparse
import colorsys
from typing import Dict, List, Tuple, Optional

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


def evenly_spaced_colors(n: int) -> List[str]:
    colors: List[str] = []
    for i in range(n):
        h = (i / max(1, n)) % 1.0
        s = 0.95
        l = 0.55
        r, g, b = colorsys.hls_to_rgb(h, l, s)
        colors.append(f"#{int(r * 255):02X}{int(g * 255):02X}{int(b * 255):02X}")
    return colors


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
    palette = evenly_spaced_colors(len(keys))
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
genes_color_dic, _gen_color = assign_colors(gene_list, color_groups)
print(f"  Gene list for coloring: {gene_list}")


# ---------------------------------------------------------------------------
# STEP 3: Match protein_ids to .gbff files by filename
# ---------------------------------------------------------------------------
print(f"Matching protein_ids to .gbff files in {args.input_dir}")
gbff_files = [f for f in os.listdir(args.input_dir) if f.endswith(".gbff")]
matched: Dict[str, str] = {}  # protein_id → gbff_path

for pid in fasta_protein_ids:
    for gbff_name in gbff_files:
        if gbff_name.endswith(f"{pid}.gbff"):
            matched[pid] = os.path.join(args.input_dir, gbff_name)
            break

not_found = [pid for pid in fasta_protein_ids if pid not in matched]
if not_found:
    print(
        f"  WARNING: {len(not_found)} protein_id(s) not found in any .gbff "
        f"filename and will be skipped: {', '.join(not_found)}"
    )
print(f"  Matched {len(matched)} / {len(fasta_protein_ids)} protein_ids")


# ---------------------------------------------------------------------------
# STEP 4: For each matched protein_id, parse the .gbff and extract context
# ---------------------------------------------------------------------------
results: List[Tuple[str, int, List[Tuple], str]] = []

for pid, gbff_path in matched.items():
    records = list(SeqIO.parse(gbff_path, "genbank"))
    if not records:
        print(f"  WARNING: empty .gbff for {pid}: {gbff_path}")
        continue

    center_locus_tag: Optional[str] = None
    matched_contig_genes: List[Tuple] = []
    taxid = "NA"
    species = "NA"
    strain = "NA"
    isolate = "NA"

    for record in records:
        contig_genes: List[Tuple[int, int, int, str, str]] = []
        found_on_contig = False
        local_taxid = "NA"
        local_species = "NA"
        local_strain = "NA"
        local_isolate = "NA"

        for feature in record.features:
            if feature.type == "source":
                db_xrefs = feature.qualifiers.get("db_xref", [])
                local_isolate = clean_metadata(
                    feature.qualifiers.get("isolate", ["NA"])[0]
                )
                local_strain = clean_metadata(
                    feature.qualifiers.get("strain", ["NA"])[0]
                )
                local_species = clean_metadata(
                    feature.qualifiers.get("organism", ["NA"])[0]
                )
                for xref in db_xrefs:
                    if xref.startswith("taxon:"):
                        local_taxid = xref.split(":")[1]
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

                if p_id == pid:
                    center_locus_tag = lt
                    found_on_contig = True

        if found_on_contig and len(contig_genes) > len(matched_contig_genes):
            matched_contig_genes = contig_genes
            taxid = local_taxid
            species = local_species
            strain = local_strain
            isolate = local_isolate

    if center_locus_tag is None:
        print(f"  WARNING: protein_id {pid} not found as CDS in {gbff_path}")
        continue
    if not matched_contig_genes:
        print(f"  WARNING: no contig genes found for {pid}")
        continue

    # -- reflect synteny --
    goi = next(
        (g for g in matched_contig_genes if g[4] == center_locus_tag), None
    )
    if goi is None:
        print(f"  WARNING: center locus_tag {center_locus_tag} not in contig genes")
        continue

    genes = list(matched_contig_genes)
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

    # -- context window (50 kb each side of the gene) --
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
if not results:
    raise SystemExit("No organisms with valid context. Nothing to write.")

label = f"Synteny_{fasta_protein_ids[0]}"
out_path = args.itol_synteny_file

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

print(f"\nWrote {len(results)} organism(s) to {out_path}")
