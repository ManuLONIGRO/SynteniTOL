#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 13:36:00 2025

UPDATED FOR NEXTFLOW
Take .gbff files from a directory and rename them based on organism and isolate info from a TSV file.
@author: mlonigro
"""
import argparse
import csv
import os
import re
import shutil
from typing import Dict, List


def sanitize(value: str) -> str:
    """Return a lowercase alphanumeric version of the given value."""
    return re.sub(r"[^a-zA-Z0-9]", "", value.strip().lower())


def load_assembly_info(tsv_path: str) -> Dict[str, Dict[str, str]]:
    """Parse TSV rows into a dict keyed by nucleotide_accession."""
    required_columns = {"nucleotide_accession", "organism"}
    info: Dict[str, Dict[str, str]] = {}

    with open(tsv_path, newline="") as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")
        missing = required_columns - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in TSV: {', '.join(sorted(missing))}")

        for row in reader:
            assembly = row["nucleotide_accession"].strip()
            if not assembly:
                continue

            info[assembly] = {
                "organism": sanitize(row["organism"]),
                "isolate": sanitize(row.get("isolate", "")),
                "protein_id": row.get("protein", ""),
            }

    return info


def resolve_gbff_inputs(path: str) -> List[str]:
    """Return a list of gbff files from either a directory or single file."""
    if os.path.isdir(path):
        gbff_files = [
            os.path.join(path, f)
            for f in os.listdir(path)
            if f.lower().endswith(".gbff")
        ]
        if not gbff_files:
            raise FileNotFoundError(f"No .gbff files found in directory {path}")
        return sorted(gbff_files)

    if os.path.isfile(path) and path.lower().endswith(".gbff"):
        return [path]

    raise FileNotFoundError(f"{path} is not a .gbff file or directory")


def build_new_name(organism: str, isolate: str, assembly: str, protein_id: str) -> str:
    """Compose the destination filename."""
    base = organism
    if isolate:
        base += isolate
    if not base:
        base = sanitize(assembly)
    return f"{base}{protein_id}.gbff".replace("_","")

def copy_with_collision_handling(src: str, dst_dir: str, new_name: str) -> str:
    """Copy src -> dst_dir/new_name, adding suffixes if needed."""
    candidate = new_name
    counter = 1
    dst_path = os.path.join(dst_dir, candidate)
    while os.path.exists(dst_path):
        candidate = f"{os.path.splitext(new_name)[0]}_{counter}.gbff"
        dst_path = os.path.join(dst_dir, candidate)
        counter += 1

    shutil.copy2(src, dst_path)
    return dst_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv", required=True, help="Path to TSV with assembly info")
    parser.add_argument("--gbff_no_assemblies", required=True, help="Path to a .gbff file (or directory containing them)")
    parser.add_argument("--output_dir_genomes", required=True, help="Output directory for renamed .gbff files")
    args = parser.parse_args()

    tsv_path = args.tsv
    gbff_input = args.gbff_no_assemblies
    output_dir = args.output_dir_genomes

    os.makedirs(output_dir, exist_ok=True)
    assembly_info = load_assembly_info(tsv_path)
    gbff_files = resolve_gbff_inputs(gbff_input)

    for gbff_file in gbff_files:
        assembly = os.path.splitext(os.path.basename(gbff_file))[0]
        info = assembly_info.get(assembly)
        if not info:
            print(f"[!] No TSV entry found for {assembly}, skipping {gbff_file}")
            continue

        new_name = build_new_name(info["organism"], info["isolate"], assembly,info["protein_id"])
        dst_file = copy_with_collision_handling(gbff_file, output_dir, new_name)
        print(f"[+] Copy {gbff_file} -> {dst_file}")


if __name__ == "__main__":
    main()
