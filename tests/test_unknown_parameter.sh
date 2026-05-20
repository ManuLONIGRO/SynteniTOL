#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TMP_LOG="$(mktemp)"
trap 'rm -f "$TMP_LOG"' EXIT

UNKNOWN_PARAM="custom_profiles"

set +e
nextflow run main.nf "--${UNKNOWN_PARAM}" profile_x.hmm >"$TMP_LOG" 2>&1
exit_code=$?
set -e

if [[ $exit_code -eq 0 ]]; then
  echo "FAIL: pipeline succeeded with unknown parameter --${UNKNOWN_PARAM}"
  cat "$TMP_LOG"
  exit 1
fi

log_content="$(<"$TMP_LOG")"
if [[ ! "$log_content" =~ Unknown\ parameter\(s\):\ --${UNKNOWN_PARAM} ]]; then
  echo "FAIL: expected unknown parameter error message not found"
  printf '%s\n' "$log_content"
  exit 1
fi

if [[ ! "$log_content" =~ Allowed\ parameters\ are: ]] || [[ ! "$log_content" =~ --custom_hmm_profiles ]]; then
  echo "FAIL: expected list of allowed parameters was not printed"
  printf '%s\n' "$log_content"
  exit 1
fi

echo "PASS: unknown CLI parameter is rejected and allowed parameters are shown"
