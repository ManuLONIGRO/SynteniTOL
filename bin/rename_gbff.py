#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 13:36:00 2025

UPDATED FOR NEXTFLOW

Copies the .gbff of every assembly to the output directory, renaming it as
{organism}{isolate}{protein_id}.gbff.

The input TSV has one row per protein of interest, so several rows (paralogs)
can point to the SAME assembly. Copying the genome once per row produced N
byte-identical copies of the same .gbff, which were then parsed N times by
the downstream steps (dataframe with repeated rows, duplicated iTOL leaves,
N identical proteomes and HMM searches).

Now the rows are grouped by assembly: the genome is copied ONCE and the
redundant copies are reported in the audit file (--out_audit).

The genome is placed in the output directory as a HARD LINK of the downloaded
file (falling back to a real copy when the link is not possible, e.g. files on
different filesystems). The bytes of every assembly are then stored only once
in work/, no matter how many times the file is needed downstream.

@author: mlonigro
"""
import os
import csv
import shutil
import argparse
import re

AUDIT_FIELDS = [
    "assembly",
    "organism",
    "isolate",
    "status",
    "n_proteins_of_interest",
    "n_copies_avoided",
    "genome_file",
    "proteins_merged",
    "copies_skipped",
]

parser = argparse.ArgumentParser()
parser.add_argument("--tsv", required=True, help="Path to input TSV file with assembly info")
parser.add_argument("--result_zip_ncbi_dir", required=True, help="Path to the result of .zip ncbi_dataset download")
parser.add_argument("--output_dir_genomes", required=True, help="Path to output directory for renamed .gbff files")
parser.add_argument("--out_audit", required=False, help="Path to the TSV audit with the deduplicated genomes")
parser.add_argument("--no_hardlink", action="store_true",
                    help="Always copy the .gbff files instead of hard-linking them (uses more disk space)")
args = parser.parse_args()

tsv_path = args.tsv
base_dir = args.result_zip_ncbi_dir
output_dir = args.output_dir_genomes
out_audit = args.out_audit
use_hardlink = not args.no_hardlink


def place_genome(src_file, dst_file):
    """Put the genome in dst_file, sharing its bytes with src_file if possible.

    A hard link keeps a single copy of the assembly bytes in work/ (the
    download directory, the renamed genomes and the merged genomes_all
    directory all point to the same inode). When the link cannot be created
    (e.g. different filesystems, or a filesystem without hard link support)
    the genome is copied as a fallback, so the pipeline always works.

    Returns "link" or "copy".
    """
    if os.path.lexists(dst_file):
        # A file with the same name already exists (e.g. two assemblies with the
        # same organism/isolate/protein). It is removed instead of overwritten:
        # overwriting would truncate it in place and, if it is a hard link to
        # another genome, it would corrupt that genome too.
        os.unlink(dst_file)
    if use_hardlink:
        # Follow symlinks: the link must point to the real file and not to the
        # symlink itself.
        real_src = os.path.realpath(src_file)
        try:
            os.link(real_src, dst_file)
            return "link"
        except OSError as exc:
            print(f"[!] Hard link not possible for {real_src} ({exc.strerror}); "
                  f"falling back to a copy.")
    shutil.copy2(src_file, dst_file)
    return "copy"


os.makedirs(output_dir, exist_ok=True)

# Group the rows by assembly: the same genome can be pointed by several
# proteins of interest (paralogs). The genome is copied only once.
rows_by_assembly = {}
with open(tsv_path, newline="") as tsvfile:
    reader = csv.DictReader(tsvfile, delimiter="\t")
    for row in reader:
        assembly = (row.get("assembly") or "").strip()
        rows_by_assembly.setdefault(assembly, []).append(row)

audit_rows = []
n_copies = 0
n_copies_avoided = 0
n_links = 0
n_real_copies = 0

for assembly, rows in rows_by_assembly.items():
    src_dir = os.path.join(base_dir, assembly)
    if not os.path.isdir(src_dir):
        # The assembly is not in this chunk: it belongs to another chunk of the
        # download (or it was never downloaded). It is reported by its own
        # chunk, so it is not added to the audit here.
        print(f"[!] The dir: {src_dir} it doesn't exists.")
        continue

    # Search for .gbff file within the directory
    gbff_files = [f for f in os.listdir(src_dir) if f.endswith(".gbff")]
    if not gbff_files:
        print(f"[!] There isn't .gbff in {src_dir}")
        audit_rows.append({
            "assembly": assembly,
            "organism": (rows[0].get("organism") or "").strip(),
            "isolate": (rows[0].get("isolate") or "").strip(),
            "status": "no_gbff",
            "n_proteins_of_interest": len(rows),
            "n_copies_avoided": 0,
            "genome_file": "",
            "proteins_merged": ";".join((r.get("protein") or "").strip() for r in rows),
            "copies_skipped": "",
        })
        continue

    src_file = os.path.join(src_dir, gbff_files[0])

    kept_file = ""
    skipped_files = []
    for row in rows:
        protein_id = (row.get("protein") or "").strip()
        organism = re.sub(r"[^a-zA-Z0-9]", "", (row.get("organism") or "").strip().lower())
        raw_isolate = row.get("isolate", "")
        isolate = re.sub(r"[^a-zA-Z0-9]", "", raw_isolate.strip().lower()) if raw_isolate.strip() else ""

        if isolate:
            new_name = f"{organism}{isolate}{protein_id}.gbff".replace("_","")
        else:
            new_name = f"{organism}{protein_id}.gbff".replace("_","")

        if not kept_file:
            # First protein of interest of this assembly: this is the copy that
            # is kept in the output directory.
            kept_file = new_name
            dst_file = os.path.join(output_dir, new_name)
            if place_genome(src_file, dst_file) == "link":
                n_links += 1
                print(f"[+] Link {src_file} -> {dst_file}")
            else:
                n_real_copies += 1
                print(f"[+] Copy {src_file} -> {dst_file}")
            n_copies += 1
        else:
            # Redundant copy of a genome already copied: it would be
            # byte-identical, so it is skipped.
            skipped_files.append(new_name)
            n_copies_avoided += 1

    if skipped_files:
        print(f"[=] {assembly}: {len(rows)} proteins of interest -> 1 genome copy "
              f"({len(skipped_files)} redundant copies avoided)")

    audit_rows.append({
        "assembly": assembly,
        "organism": (rows[0].get("organism") or "").strip(),
        "isolate": (rows[0].get("isolate") or "").strip(),
        "status": "copied",
        "n_proteins_of_interest": len(rows),
        "n_copies_avoided": len(skipped_files),
        "genome_file": kept_file,
        "proteins_merged": ";".join((r.get("protein") or "").strip() for r in rows),
        "copies_skipped": ";".join(skipped_files),
    })

if out_audit:
    out_dir = os.path.dirname(out_audit)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_audit, "w", newline="") as auditfile:
        writer = csv.DictWriter(auditfile, fieldnames=AUDIT_FIELDS, delimiter="\t")
        writer.writeheader()
        for audit_row in audit_rows:
            writer.writerow(audit_row)
    print(f"[+] Audit: {out_audit}")

print(f"[=] Genomes copied: {n_copies} | redundant copies avoided: {n_copies_avoided} "
      f"| hard links: {n_links} | real copies: {n_real_copies}")
