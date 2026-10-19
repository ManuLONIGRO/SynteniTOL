#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-genome locus_tag mapping.

Parses ONE .gbff a single time and joins it against the hmmsearch .result
files of ALL profiles for that genome. Replaces the old per-profile
invocation of ids_locustag_mapping.py, which re-parsed the same .gbff once
per profile (O(profiles x genomes)).

Output .map / .no_map files are byte-compatible with the old script:
  mapped   : "<protein_id> -> <locus_tag>\tcog=<cog>\tlen=<L>\taligned=<A>\tevalue=<E>\thmmsearch_acc=<X>\n"
  no hits  : "no hits\n"
  no_map   : "<protein_id>\n"
  gbff missing / result missing or unparsable : write nothing (parity)

@author: mlonigro
"""

import os
import sys
import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--genomes_dir", required=True, help="Path to the directory containing .gbff files")
parser.add_argument("--genome", required=True, help="Genome/proteome base name (without .gbff/.faa)")
parser.add_argument("--result_files", nargs="+", required=True, help="hmmsearch .result files of this genome (one per profile)")
parser.add_argument("--gene_names", required=True, help="Comma-separated profile labels aligned with --result_files")
parser.add_argument("--out_dir", required=True, help="Path to the output directory for .map/.no_map files")
parser.add_argument("--evalue", type=float, default=1e-10, help="E-value threshold for filtering")
args = parser.parse_args()


def parse_result(result_path, evalue_limit):
    """Parse one hmmsearch domtblout file.

    Returns:
      None  -> file missing or unparsable (write nothing, parity with old script)
      {}    -> no hit passes the e-value filter (write "no hits")
      dict  -> protein_id -> {total_aligned, evalue, hmm_acc}
    """
    hit_stats = {}
    try:
        with open(result_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.split()
                if len(parts) < 19:
                    # strange/incomplete line; skip it
                    continue

                # full-sequence E-value is column 7 (1-based) -> parts[6]
                evalue = float(parts[6])
                if evalue > evalue_limit:
                    continue

                hit_id = parts[0]
                ali_from = int(parts[17])
                ali_to = int(parts[18])
                domain_len = ali_to - ali_from + 1

                stats = hit_stats.setdefault(
                    hit_id,
                    {"total_aligned": 0, "evalue": evalue, "hmm_acc": float(parts[21])},
                )
                stats["total_aligned"] += domain_len
                # Old script overwrites hmm_acc on every line => last value wins.
                stats["hmm_acc"] = float(parts[21])
                # Save the best (minimum) e-value observed for that hit
                if evalue < stats["evalue"]:
                    stats["evalue"] = evalue
    except FileNotFoundError:
        return None
    except Exception as e:  # noqa: BLE001 - mirror old behaviour (skip this file)
        print(f"Error parsing {result_path}: {e}", file=sys.stderr)
        return None

    return hit_stats


def main():
    cogs = args.gene_names.split(",")
    if len(cogs) != len(args.result_files):
        print(
            f"ERROR: {len(cogs)} gene names vs {len(args.result_files)} result files",
            file=sys.stderr,
        )
        sys.exit(1)

    # 1) The .gbff is always checked first: if missing, write nothing at all
    #    (same behaviour as the old per-profile script).
    gbff_file = os.path.join(args.genomes_dir, args.genome + ".gbff")
    if not os.path.exists(gbff_file):
        print(f"GBFF not found for {args.genome}: {gbff_file}", file=sys.stderr)
        return

    # 2) Parse every .result file of this genome (cheap) and collect the
    #    protein_ids that passed the filter.
    entries = []  # (cog, result_path, stats)
    needed = set()
    for cog, result_path in zip(cogs, args.result_files):
        stats = parse_result(result_path, args.evalue)
        entries.append((cog, result_path, stats))
        if stats:
            needed.update(stats.keys())

    # 3) Parse the .gbff ONCE (only if at least one hit needs mapping).
    #    protein_id -> [(locus_tag, length), ...]   (list preserves duplicates)
    pid2info = {}
    if needed:
        try:
            for record in SeqIO.parse(gbff_file, "genbank"):
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    qualifiers = feature.qualifiers
                    protein_id = qualifiers.get("protein_id", [""])[0]
                    if not protein_id or protein_id not in needed:
                        continue
                    locus_tag = qualifiers.get("locus_tag", [""])[0] or protein_id
                    length_seq = len(qualifiers.get("translation", [""])[0])
                    pid2info.setdefault(protein_id, []).append((locus_tag, length_seq))
        except Exception as e:  # noqa: BLE001
            print(f"Error parsing {gbff_file}: {e}", file=sys.stderr)

    # 4) Emit one .map / .no_map per profile, identical format to old script.
    os.makedirs(args.out_dir, exist_ok=True)
    for cog, result_path, stats in entries:
        base = os.path.basename(result_path)
        if base.endswith(".result"):
            base = base[: -len(".result")]
        map_path = os.path.join(args.out_dir, base + ".map")
        no_map_path = os.path.join(args.out_dir, base + ".no_map")

        if stats is None:  # error/missing -> write nothing (parity)
            continue

        # Old script always opens both files (creating an empty .no_map) as
        # long as there is something to report. Replicate that.
        with open(map_path, "w") as out_map, open(no_map_path, "w") as out_no_map:
            if not stats:
                out_map.write("no hits\n")
                continue

            for hit_id, s in stats.items():
                infos = pid2info.get(hit_id)
                if infos:
                    for locus_tag, length_seq in infos:
                        out_map.write(
                            f"{hit_id} -> {locus_tag}\tcog={cog}\tlen={length_seq}\t"
                            f"aligned={s['total_aligned']}\tevalue={s['evalue']}\t"
                            f"hmmsearch_acc={s['hmm_acc']}\n"
                        )
                else:
                    out_no_map.write(f"{hit_id}\n")


if __name__ == "__main__":
    main()
