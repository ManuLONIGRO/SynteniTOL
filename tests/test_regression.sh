#!/usr/bin/env bash
set -euo pipefail

# Regression test for SynteniToL.
#
# 1. Replays the exact command stored in tests/test_result/test_run_command.txt
#    (produced by the stable run results/run_20260916_142328) against the
#    current version of the pipeline, using tests/test_result as golden set.
# 2. Compares the newly produced itol_*.txt files against the golden reference
#    with tests/compare_itol_results.py.
# 3. Prints a WARNING (and a detailed diff report) if anything changed, so the
#    user can review whether the differences are expected after the last
#    implementation.
#
# The fresh run is kept on disk (default: results/regression_<timestamp>) so it
# can be inspected and re-compared afterwards.
#
# Usage:
#   tests/test_regression.sh [--outdir <path>]
#   tests/test_regression.sh --help
# Any other option aborts before the pipeline runs (exit 2).
#
# Exit codes:
#   0 = results identical to the golden reference
#   1 = pipeline or comparison completed but differences found (WARNING)
#   2 = pipeline failed, comparison errored, or an unknown option was given

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

GOLDEN_DIR="tests/test_result"
COMMAND_FILE="${GOLDEN_DIR}/test_run_command.txt"

# Resolve to an absolute path (also works if GOLDEN_DIR is given as absolute).
if [[ "$GOLDEN_DIR" == /* ]]; then
    reference_dir="$GOLDEN_DIR"
else
    reference_dir="$ROOT_DIR/$GOLDEN_DIR"
fi
COMMAND_FILE="$reference_dir/test_run_command.txt"

USAGE="Uso: tests/test_regression.sh [--outdir <path>]"
AVAILABLE_OPTS="Opciones disponibles:
  --outdir <path>   directorio de salida del nuevo run (default: results/regression_<timestamp>)
  -h, --help        mostrar este mensaje y salir"

OUTDIR="results/regression_$(date +%Y%m%d_%H%M%S)"

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

if [[ ! -f "$COMMAND_FILE" ]]; then
    echo "FAIL: golden reference file not found: $COMMAND_FILE" >&2
    exit 2
fi

command_str="$(<"$COMMAND_FILE")"
command_str="${command_str} --outdir ${OUTDIR}"

echo "=== Test de regresion de SynteniToL ==="
echo "  Golden de referencia : ${reference_dir}/"
echo "  Run nuevo            : ${OUTDIR}/"
echo "  Comando a ejecutar   : ${command_str}"
echo
echo "NOTA: el test re-descarga genomas desde NCBI y puede tardar ~10 min."
echo "      Si NCBI actualizo un ensamblaje, el diff lo mostrara como cambio."
echo

mkdir -p "$(dirname "$OUTDIR")"

# Export so a custom test command (or a simulation) can reference the outdir.
export OUTDIR

TMP_LOG="$(mktemp)"
trap 'rm -f "$TMP_LOG"' EXIT

set +e
bash -c "$command_str" >"$TMP_LOG" 2>&1
pipeline_exit=$?
set -e

if [[ $pipeline_exit -ne 0 ]]; then
    echo "FAIL: el pipeline no termino correctamente (exit code ${pipeline_exit})." >&2
    echo "Ultimas lineas del log de Nextflow:" >&2
    tail -50 "$TMP_LOG" >&2
    rm -rf "$OUTDIR"
    exit 2
fi

set +e
python3 "$ROOT_DIR/tests/compare_itol_results.py" \
    --reference "$reference_dir" \
    --current "$OUTDIR" \
    --diff-out "${OUTDIR}/itol_diff_report.txt"
compare_exit=$?
set -e

case $compare_exit in
    0)
        echo
        echo "PASS: los 6 archivos itol* del nuevo run coinciden con la referencia."
        exit 0
        ;;
    1)
        echo
        echo "WARNING: el pipeline produce resultados DIFERENTES a la referencia."
        echo "Revisa el diff detallado en: ${OUTDIR}/itol_diff_report.txt"
        echo "y decide si los cambios corresponden a la ultima implementacion."
        exit 1
        ;;
    *)
        echo
        echo "WARNING: no se pudieron comparar los resultados (error de comparacion)." >&2
        echo "Revisa manualmente con:" >&2
        echo "  python3 ${ROOT_DIR}/tests/compare_itol_results.py --reference ${reference_dir} --current ${OUTDIR}" >&2
        exit 2
        ;;
esac