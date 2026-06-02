#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TMP_LOG="$(mktemp)"
TMP_FASTA="$(mktemp)"
TMP_DIR="$(mktemp -d)"
trap 'rm -f "$TMP_LOG" "$TMP_FASTA"; rm -rf "$TMP_DIR"' EXIT

printf '>dummy\nACGT\n' >"$TMP_FASTA"
mkdir -p "${TMP_DIR}/subdir"
touch "${TMP_DIR}/profile_a.hmm" "${TMP_DIR}/subdir/profile_a.hmm"

run_pipeline() {
  local custom_profiles="$1"
  local expected_pattern="$2"

  : >"$TMP_LOG"
  set +e
  nextflow run main.nf \
    --inputFASTA "$TMP_FASTA" \
    --cogs COG1152,COG1795 \
    --custom_hmm_profiles "$custom_profiles" \
    >"$TMP_LOG" 2>&1
  local exit_code=$?
  set -e

  if [[ $exit_code -eq 0 ]]; then
    echo "FAIL: pipeline succeeded with custom profiles: ${custom_profiles}"
    cat "$TMP_LOG"
    exit 1
  fi

  local log_content
  log_content="$(<"$TMP_LOG")"
  if [[ ! "$log_content" =~ ${expected_pattern} ]]; then
    echo "FAIL: expected error matching '${expected_pattern}' not found"
    printf '%s\n' "$log_content"
    exit 1
  fi
}

run_pipeline \
  "${TMP_DIR}/profile_a.hmm,${TMP_DIR}/profile_a.hmm" \
  'Duplicate paths in --custom_hmm_profiles:'

run_pipeline \
  "${TMP_DIR}/profile_a.hmm,${TMP_DIR}/subdir/profile_a.hmm" \
  'Duplicate custom profile names in --custom_hmm_profiles:'

echo "PASS: duplicate custom HMM profiles are rejected before the pipeline runs"
