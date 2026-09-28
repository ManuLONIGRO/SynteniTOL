#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TMP_LOG="$(mktemp)"
TMP_FASTA="$(mktemp)"
TMP_OUT="$(mktemp -d)"
trap 'rm -f "$TMP_LOG" "$TMP_FASTA"; rm -rf "$TMP_OUT"' EXIT

: > "$TMP_FASTA"

# --- Case A: unknown gene in --gene_order must stop the pipeline ---
set +e
nextflow run main.nf \
  --inputFASTA "$TMP_FASTA" \
  --cogs COG1152,COG1795 \
  --gene_order COG1152,COG9999 \
  --outdir "$TMP_OUT/unknown_gene" \
  >"$TMP_LOG" 2>&1
exit_code=$?
set -e

if [[ $exit_code -eq 0 ]]; then
  echo "FAIL: pipeline succeeded with an unknown gene in --gene_order"
  cat "$TMP_LOG"
  exit 1
fi

log_content="$(<"$TMP_LOG")"
if [[ ! "$log_content" =~ Invalid\ identifiers\ in\ --gene_order: ]]; then
  echo "FAIL: expected 'Invalid identifiers in --gene_order' error message not found"
  printf '%s\n' "$log_content"
  exit 1
fi

if [[ ! "$log_content" =~ COG9999 ]]; then
  echo "FAIL: expected unknown gene COG9999 to be reported"
  printf '%s\n' "$log_content"
  exit 1
fi

echo "PASS: unknown genes in --gene_order are rejected before the pipeline runs"

# --- Case B: duplicates in --gene_order must stop the pipeline ---
set +e
nextflow run main.nf \
  --inputFASTA "$TMP_FASTA" \
  --cogs COG1152,COG1795 \
  --gene_order COG1795,COG1152,COG1795 \
  --outdir "$TMP_OUT/duplicate_genes" \
  >"$TMP_LOG" 2>&1
exit_code=$?
set -e

if [[ $exit_code -eq 0 ]]; then
  echo "FAIL: pipeline succeeded with duplicate genes in --gene_order"
  cat "$TMP_LOG"
  exit 1
fi

log_content="$(<"$TMP_LOG")"
if [[ ! "$log_content" =~ Duplicate\ identifiers\ in\ --gene_order: ]]; then
  echo "FAIL: expected duplicate --gene_order error message not found"
  printf '%s\n' "$log_content"
  exit 1
fi

echo "PASS: duplicates in --gene_order are rejected before the pipeline runs"

# --- Case C: genes missing from --gene_order only warn, run continues ---
# An empty FASTA is used so the warning is emitted (it happens before any
# process) and the pipeline then fails deterministically on empty assemblies
# (.ifEmpty in main.nf) without contacting NCBI. The timeout -k guarantees the
# test never hangs even if something unexpected blocks.
set +e
timeout -k 60 180 nextflow run main.nf \
  --inputFASTA "$TMP_FASTA" \
  --cogs COG1152,COG1795 \
  --gene_order COG1795 \
  --outdir "$TMP_OUT/partial_list" \
  >"$TMP_LOG" 2>&1
exit_code=$?
set -e

log_content="$(<"$TMP_LOG")"
if [[ "$log_content" =~ Invalid\ identifiers\ in\ --gene_order: ]]; then
  echo "FAIL: valid partial --gene_order was rejected as invalid"
  printf '%s\n' "$log_content"
  exit 1
fi

if [[ ! "$log_content" =~ missing\ in\ --gene_order: ]]; then
  echo "FAIL: expected 'missing in --gene_order' warning message not found"
  printf '%s\n' "$log_content"
  exit 1
fi

if [[ ! "$log_content" =~ COG1152 ]]; then
  echo "FAIL: expected missing gene COG1152 to be reported in the warning"
  printf '%s\n' "$log_content"
  exit 1
fi

echo "PASS: genes missing from --gene_order produce a warning without stopping the pipeline"
