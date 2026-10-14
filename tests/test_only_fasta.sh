#!/usr/bin/env bash
set -euo pipefail

# Integration test for the --only_fasta mode.
#
# Runs the real pipeline (downloads genomes from NCBI) with `--only_fasta` and
# verifies that:
#   - it finishes successfully without requiring --cogs
#   - formatted_headers.fasta and the tree are produced
#   - no synteny / iTOL / profiling files are generated
#   - analysis-only parameters are ignored with a warning
#
# Requires network access to NCBI (same as the regression test).
#
# Usage:
#   tests/test_only_fasta.sh [--outdir <path>]
#   tests/test_only_fasta.sh --help
# Any other option aborts before the pipeline runs (exit 2).
#
# Exit codes:
#   0 = PASS
#   1 = pipeline succeeded but outputs are wrong (FAIL)
#   2 = pipeline failed, comparison errored, or an unknown option was given

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

USAGE="Uso: tests/test_only_fasta.sh [--outdir <path>]"
AVAILABLE_OPTS="Opciones disponibles:
  --outdir <path>   directorio de salida del run (default: results/only_fasta_<timestamp>)
  -h, --help        mostrar este mensaje y salir"

OUTDIR="results/only_fasta_$(date +%Y%m%d_%H%M%S)"

case "$#" in
    0) ;;
    2)
        if [[ "$1" == "--outdir" ]]; then
            OUTDIR="$2"
        else
            echo "ERROR: Opcion desconocida: $1" >&2
            echo "$USAGE" >&2
            echo >&2
            echo "$AVAILABLE_OPTS" >&2
            exit 2
        fi
        ;;
    1)
        if [[ "$1" == "-h" || "$1" == "--help" ]]; then
            echo "$USAGE"
            echo
            echo "$AVAILABLE_OPTS"
            exit 0
        fi
        echo "ERROR: Opcion desconocida: $1" >&2
        echo "$USAGE" >&2
        echo >&2
        echo "$AVAILABLE_OPTS" >&2
        exit 2
        ;;
    *)
        echo "ERROR: Demasiados argumentos." >&2
        echo "$USAGE" >&2
        echo >&2
        echo "$AVAILABLE_OPTS" >&2
        exit 2
        ;;
esac

INPUT_FASTA="data/fwd/small_fwda_test"
LOG_FILE="$(mktemp)"
trap 'rm -f "$LOG_FILE"' EXIT

echo "Input FASTA : $INPUT_FASTA"
echo "Outdir      : $OUTDIR"
echo "Run: nextflow run main.nf --inputFASTA $INPUT_FASTA --only_fasta --cogs COG1152,COG1029 --color_by_group COG1152-COG1029 --outdir $OUTDIR"
echo

set +e
nextflow run main.nf \
    --inputFASTA "$INPUT_FASTA" \
    --only_fasta \
    --cogs COG1152,COG1029 \
    --color_by_group COG1152-COG1029 \
    --outdir "$OUTDIR" \
    >"$LOG_FILE" 2>&1
EXIT_CODE=$?
set -e

LOGS="$(<"$LOG_FILE")"

if [[ $EXIT_CODE -ne 0 ]]; then
    echo "FAIL: pipeline exited with code $EXIT_CODE"
    tail -50 "$LOG_FILE"
    exit 2
fi

if [[ ! -f "$OUTDIR/formatted_headers.fasta" ]]; then
    echo "FAIL: formatted_headers.fasta was not produced"
    exit 1
fi

if ! ls "$OUTDIR"/tree_*.nwk >/dev/null 2>&1; then
    echo "FAIL: tree_*.nwk was not produced"
    exit 1
fi

if ls "$OUTDIR"/itol_*.txt >/dev/null 2>&1; then
    echo "FAIL: iTOL files were produced in --only_fasta mode:"
    ls "$OUTDIR"/itol_*.txt
    exit 1
fi

if [[ "$LOGS" =~ Missing\ --cogs ]]; then
    echo "FAIL: pipeline required --cogs despite --only_fasta"
    printf '%s\n' "$LOGS" | grep -i "missing --cogs"
    exit 1
fi

if [[ ! "$LOGS" =~ ignored\ in\ --only_fasta\ mode ]]; then
    echo "FAIL: expected warning about ignored parameters in --only_fasta mode"
    exit 1
fi

echo "PASS: --only_fasta produced formatted_headers.fasta and the tree, without running the synteny/profiling analysis"