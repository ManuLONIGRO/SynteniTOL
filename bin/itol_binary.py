#!/usr/bin/env python3
"""
Build organism-level gene presence from the synteny dataframe.

Groups hits by organism (locus_tag prefix before '_'), maps each organism to the
same header format used by formatted_headers.fasta and itol_synteny (via
df_complete.tsv), and writes one row per organism with its gene list.
"""

import argparse

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--dataframe", type=str, required=True, help="Path to the input dataframe TSV file")
parser.add_argument("--df_complete", type=str, required=True, help="Path to df_complete.tsv with organism metadata")
parser.add_argument("--out_tsv", type=str, required=True, help="Path to the output TSV file (organism, genes)")
parser.add_argument("--acc_threshold", type=float, default=0.8, help="Minimum accuracy (acc) to keep a hit")
args = parser.parse_args()


def sanitize(value) -> str:
    if pd.isna(value):
        return "NA"
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return "NA"
    return (
        text
        .replace(" ", "_")
        .replace(";", "_")
        .replace(":", "_")
        .replace("=", "_")
        .replace("(", "_")
        .replace(")", "_")
    )


def build_headers_from_df_complete(df_complete: pd.DataFrame) -> tuple[dict[str, str], dict[str, str]]:
    """Return protein_id -> header and locus prefix -> header mappings."""
    header_count: dict[str, int] = {}
    protein_to_header: dict[str, str] = {}
    prefix_to_header: dict[str, str] = {}

    for _, row in df_complete.iterrows():
        taxid = sanitize(row.get("local_taxid", "NA"))
        species = sanitize(row.get("local_species", "NA"))
        strain = sanitize(row.get("local_strain", "NA"))
        isolate = sanitize(row.get("local_isolate", "NA"))

        base = f"{taxid}|{species}|{strain}|{isolate}"
        header_count[base] = header_count.get(base, 0) + 1
        header = f"{base}|{header_count[base]}"

        protein_id = str(row.get("protein_id", "")).strip()
        if protein_id and protein_id.lower() != "nan":
            protein_to_header[protein_id] = header

        locus_tag = row.get("locus_tag")
        if pd.notna(locus_tag):
            prefix = str(locus_tag).split("_")[0]
            prefix_to_header[prefix] = header

    return protein_to_header, prefix_to_header


def resolve_organism_name(
    organism: str,
    protein_ids: list[str],
    protein_to_header: dict[str, str],
    prefix_to_header: dict[str, str],
) -> str:
    for protein_id in protein_ids:
        if protein_id in protein_to_header:
            return protein_to_header[protein_id]
    return prefix_to_header.get(organism, organism)


print(f"Reading dataframe from {args.dataframe}")
df = pd.read_csv(
    args.dataframe,
    sep="\t",
    dtype={"locus_tag": str, "gene": str, "protein_id": str},
)
df["acc"] = pd.to_numeric(df["acc"], errors="coerce")
df_filtered = df[df["acc"] >= args.acc_threshold].copy()
df_filtered["organism"] = df_filtered["locus_tag"].str.split("_").str[0]

print(f"Reading organism metadata from {args.df_complete}")
df_complete = pd.read_csv(
    args.df_complete,
    sep="\t",
    dtype={"protein_id": str, "locus_tag": str, "local_taxid": str},
)
protein_to_header, prefix_to_header = build_headers_from_df_complete(df_complete)

organism_genes = (
    df_filtered.groupby("organism")
    .agg(
        genes=("gene", lambda values: sorted(set(values))),
        protein_ids=("protein_id", lambda values: list(dict.fromkeys(str(v) for v in values))),
    )
    .reset_index()
)
organism_genes["organism_name"] = organism_genes.apply(
    lambda row: resolve_organism_name(
        row["organism"],
        row["protein_ids"],
        protein_to_header,
        prefix_to_header,
    ),
    axis=1,
)

print(f"Writing organism gene presence for {len(organism_genes)} organisms to {args.out_tsv}")
with open(args.out_tsv, "w") as out:
    out.write("organism\tgenes\n")
    for _, row in organism_genes.iterrows():
        genes_str = ",".join(row["genes"])
        out.write(f"{row['organism_name']}\t{genes_str}\n")
