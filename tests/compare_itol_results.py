#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare the iTOL output files (itol_*.txt) of two directories: a golden/reference
directory and a freshly produced run. It reports every difference so the user can
review whether the changes make sense after an implementation.

Usage:
  python3 tests/compare_itol_results.py \
      --reference tests/test_result \
      --current results/run_YYYYMMDD_HHMMSS \
      [--diff-out /path/to/itol_diff_report.txt]

Exit codes:
  0 = both directories produce the same itol* files
  1 = differences found (warning); a detailed diff report was written
  2 = usage error, missing reference/current dir, or a referenced file is missing

Notes:
  - Golden files usually carry a test_ prefix (e.g. test_itol_profiling.txt) so
    they are not affected by the *itol_*.txt* entry in .gitignore. The prefix is
    detected automatically, so this script also works on two raw results dirs.
  - Data rows are compared in sorted order so Nextflow parallelism does not
    produce spurious changes. A change that only reorders rows is reported as
    "order only".
"""

import argparse
import datetime
import difflib
import glob
import os
import sys
from collections import Counter, defaultdict

# Header lines that are purely cosmetic (datasets labels / legend titles).
COSMETIC_HEADER = {"DATASET_LABEL", "LEGEND_TITLE"}


def eprint(msg):
    print(msg, file=sys.stderr)


def error(msg):
    eprint(f"ERROR: {msg}")
    sys.exit(2)


def read_lines(path):
    with open(path, "r", encoding="utf-8") as fh:
        return [line.rstrip("\r\n") for line in fh]


def split_header_data(lines):
    """Split file lines into a header block (up to and including DATA) and the
    remaining data rows."""
    try:
        idx = lines.index("DATA")
    except ValueError:
        return lines, []
    return lines[: idx + 1], lines[idx + 1 :]


def find_itol_files(directory):
    """Return all itol_*.txt (or test_itol_*.txt) files in a directory."""
    files = sorted(glob.glob(os.path.join(directory, "test_itol_*.txt")))
    if not files:
        files = sorted(glob.glob(os.path.join(directory, "itol_*.txt")))
    return files


def relative_name(reference_path):
    """Map test_itol_X.txt in the golden dir to itol_X.txt in the run dir."""
    base = os.path.basename(reference_path)
    if base.startswith("test_"):
        return base[len("test_") :]
    return base


def classify_header_line(line):
    key = line.split(",", 1)[0]
    return "cosmetico" if key in COSMETIC_HEADER else "estructural"


def organism_key(row):
    return row.split("|", 1)[0]


def row_summary(ref_rows, cur_rows):
    """Return a short, human readable summary of which organisms changed."""
    ref_by_org = defaultdict(list)
    cur_by_org = defaultdict(list)
    for row in ref_rows:
        ref_by_org[organism_key(row)].append(row)
    for row in cur_rows:
        cur_by_org[organism_key(row)].append(row)

    common = sorted(set(ref_by_org) & set(cur_by_org), key=lambda o: o)
    removed_orgs = sorted(set(ref_by_org) - set(cur_by_org), key=lambda o: o)
    added_orgs = sorted(set(cur_by_org) - set(ref_by_org), key=lambda o: o)
    modified = [
        org
        for org in common
        if Counter(ref_by_org[org]) != Counter(cur_by_org[org])
    ]

    lines = []
    if removed_orgs:
        lines.append(f"    organismos solo en referencia ({len(removed_orgs)}): {removed_orgs}")
    if added_orgs:
        lines.append(f"    organismos solo en run        ({len(added_orgs)}): {added_orgs}")
    if modified:
        lines.append(f"    organismos con filas cambiadas ({len(modified)}): {modified}")
    if removed_orgs or added_orgs or modified:
        return "\n".join(lines)
    return None


def compare_file(reference_path, current_path, diff_fh):
    """Compare a single pair of files. Returns (status, detail_lines)."""
    ref_lines = read_lines(reference_path)
    cur_lines = read_lines(current_path)

    ref_header, ref_data = split_header_data(ref_lines)
    cur_header, cur_data = split_header_data(cur_lines)

    # ---- Header comparison ----
    header_issues = []
    if ref_header != cur_header:
        for diff_line in difflib.ndiff(ref_header, cur_header):
            marker = diff_line[:2]
            if marker in ("- ", "+ "):
                kind = classify_header_line(diff_line[2:])
                header_issues.append(f"    [{kind}] {diff_line[0]} {diff_line[2:]}")

    # ---- Data comparison (order normalized) ----
    ref_sorted = sorted(ref_data)
    cur_sorted = sorted(cur_data)

    data_issues = []
    order_only = False
    if ref_sorted != cur_sorted:
        ref_multiset = Counter(ref_data)
        cur_multiset = Counter(cur_data)
        if ref_multiset == cur_multiset:
            order_only = True
        summary = row_summary(ref_data, cur_data)
        if summary:
            data_issues.append(summary)

    # ---- Write detailed diff ----
    diff_fh.write(f"\n{'=' * 74}\n")
    diff_fh.write(f"Archivo: {os.path.basename(reference_path)}\n")
    diff_fh.write(f"  referencia: {reference_path}\n")
    diff_fh.write(f"  run      : {current_path}\n")
    diff_fh.write("=" * 74 + "\n")

    if header_issues:
        diff_fh.write("\nCabecera (diferencias):\n")
        for issue in header_issues:
            diff_fh.write(issue + "\n")
    if order_only:
        diff_fh.write("\nSolo cambio el ORDEN de las filas de datos (mismos contenidos).\n")
    if data_issues:
        diff_fh.write("\nDatos (resumen):\n")
        for issue in data_issues:
            diff_fh.write(issue + "\n")
    if ref_sorted != cur_sorted:
        diff_fh.write("\nDiff de filas de datos (ordenadas):\n")
        for line in difflib.unified_diff(
            ref_sorted,
            cur_sorted,
            fromfile=os.path.join("referencia", os.path.basename(reference_path)),
            tofile=os.path.join("run", os.path.basename(current_path)),
            lineterm="",
        ):
            diff_fh.write(line + "\n")

    if header_issues or ref_sorted != cur_sorted:
        return "CAMBIADO", header_issues, data_issues, order_only
    return "IDENTICO", header_issues, data_issues, order_only


def main():
    parser = argparse.ArgumentParser(
        description="Compare itol_*.txt outputs of a reference dir vs a new run."
    )
    parser.add_argument("--reference", required=True, help="Golden/reference directory")
    parser.add_argument("--current", required=True, help="Directory with the new run results")
    parser.add_argument(
        "--diff-out",
        default=None,
        help="Path for the detailed diff report (default: <current>/itol_diff_report.txt)",
    )
    args = parser.parse_args()

    if not os.path.isdir(args.reference):
        error(f"reference directory not found: {args.reference}")
    if not os.path.isdir(args.current):
        error(f"current (results) directory not found: {args.current}")

    ref_files = find_itol_files(args.reference)
    if not ref_files:
        error(f"no itol_*.txt files found in reference directory: {args.reference}")

    diff_out = args.diff_out or os.path.join(args.current, "itol_diff_report.txt")
    with open(diff_out, "w", encoding="utf-8") as diff_fh:
        diff_fh.write("Comparacion de resultados itol* (referencia vs run)\n")
        diff_fh.write(f"  referencia : {args.reference}\n")
        diff_fh.write(f"  run        : {args.current}\n")
        diff_fh.write(f"  fecha      : {datetime.datetime.now().isoformat()}\n")

        n_identical = 0
        problems = []  # (name, description, is_error)
        for ref_file in ref_files:
            name = relative_name(ref_file)
            cur_file = os.path.join(args.current, name)
            if not os.path.isfile(cur_file):
                cur_file = os.path.join(args.current, "test_" + name)
            if not os.path.isfile(cur_file):
                problems.append((name, "archivo esperado falta en el run", True))
                continue
            status, header_issues, data_issues, order_only = compare_file(
                ref_file, cur_file, diff_fh
            )
            if status == "IDENTICO":
                n_identical += 1
                print(f"[OK]      {name}: IDENTICO")
            else:
                parts = []
                if header_issues:
                    parts.append(f"{len(header_issues)} lineas de cabecera")
                if order_only:
                    parts.append("solo orden de filas")
                elif data_issues:
                    parts.append("filas de datos")
                problems.append((name, ", ".join(parts) if parts else "contenido", False))
                print(f"[WARNING] {name}: CAMBIADO ({', '.join(parts)})")

        # Detect files produced by the run but absent in the reference.
        cur_itol = sorted(glob.glob(os.path.join(args.current, "itol_*.txt")))
        cur_names = {
            os.path.basename(f)
            for f in cur_itol
            if os.path.basename(f) != os.path.basename(diff_out)
        }
        ref_names = {os.path.basename(relative_name(f)) for f in ref_files}
        for extra in sorted(cur_names - ref_names):
            problems.append((extra, "archivo nuevo no presente en la referencia", False))
            print(f"[WARNING] {extra}: NUEVO (no existe en la referencia)")

        print()
    # end with block

    diff_fh_hint = f"\nDiff detallado: {diff_out}"
    if not problems and n_identical == len(ref_files):
        print(f"PASS: los {n_identical} archivos itol* coinciden con la referencia.")
        print(diff_fh_hint)
        sys.exit(0)

    print()
    print("WARNING: se detectaron diferencias frente a la referencia:")
    has_error = False
    for name, desc, is_error in problems:
        marker = "ERROR" if is_error else "diff"
        print(f"  [{marker}] {name}: {desc}")
        has_error = has_error or is_error
    print(diff_fh_hint)
    sys.exit(2 if has_error else 1)


if __name__ == "__main__":
    main()