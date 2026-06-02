#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TMP_LOG="$(mktemp)"
TMP_FASTA="$(mktemp)"
trap 'rm -f "$TMP_LOG" "$TMP_FASTA"' EXIT

printf '>dummy\nACGT\n' >"$TMP_FASTA"

set +e
nextflow run main.nf \
  --inputFASTA "$TMP_FASTA" \
  --cogs COG1152,COG1795,COG1152 \
  >"$TMP_LOG" 2>&1
exit_code=$?
set -e

if [[ $exit_code -eq 0 ]]; then
  echo "FAIL: pipeline succeeded with duplicate COGs in --cogs"
  cat "$TMP_LOG"
  exit 1
fi

log_content="$(<"$TMP_LOG")"
if [[ ! "$log_content" =~ Duplicate\ COG\ identifiers\ in\ --cogs: ]]; then
  echo "FAIL: expected duplicate COG error message not found"
  printf '%s\n' "$log_content"
  exit 1
fi

if [[ ! "$log_content" =~ COG1152 ]]; then
  echo "FAIL: expected duplicate COG1152 to be reported"
  printf '%s\n' "$log_content"
  exit 1
fi

echo "PASS: duplicate COGs in --cogs are rejected before the pipeline runs"
