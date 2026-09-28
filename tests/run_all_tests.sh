#!/usr/bin/env bash
set -uo pipefail

# Runs the full SynteniToL test suite:
#   - test_unknown_parameter.sh        (CLI validation)
#   - test_duplicate_cogs.sh           (CLI validation)
#   - test_duplicate_custom_profiles.sh (CLI validation)
#   - test_gene_order.sh               (CLI validation, ~1-3 min)
#   - test_regression.sh               (golden reference comparison, FULL pipeline ~11 min)
#
# Expected total duration: ~13-14 min (validation + regression). ~11 min of that
# is the regression re-running the real pipeline, which is normal, not a hang.
#
# Usage:
#   tests/run_all_tests.sh [--skip-regression]
#   tests/run_all_tests.sh --help
# Any other option aborts the suite (exit 2) without running any test.
#
# Exit codes:
#   0 = all tests passed
#   1 = at least one test passed with a WARNING (e.g. regression found differences)
#   2 = at least one test FAILED, a test timed out, or an unknown option was given

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

USAGE="Uso: tests/run_all_tests.sh [--skip-regression]"
AVAILABLE_OPTS="Opciones disponibles:
  --skip-regression   omitir test_regression.sh (el pipeline tarda ~11 min; total ~13-14 min)
  -h, --help          mostrar este mensaje y salir"

SKIP_REGRESSION=0
for arg in "$@"; do
    case "$arg" in
        --skip-regression)
            SKIP_REGRESSION=1
            ;;
        -h|--help)
            echo "$USAGE"
            echo
            echo "$AVAILABLE_OPTS"
            exit 0
            ;;
        *)
            echo "ERROR: Opcion desconocida: $arg" >&2
            echo "$USAGE" >&2
            echo >&2
            echo "$AVAILABLE_OPTS" >&2
            exit 2
            ;;
    esac
done

VALIDATION_TESTS=(
    test_unknown_parameter.sh
    test_duplicate_cogs.sh
    test_duplicate_custom_profiles.sh
    test_gene_order.sh
)

declare -a RESULTS=()
OVERALL=0
FAILED=0

run_test() {
    local name="$1"
    set +e
    timeout -k 60 600 bash "tests/${name}" >/tmp/syntenitol_test_${name}.log 2>&1
    local code=$?
    set -e
    RESULTS+=("${name}:${code}")
    return $code
}

echo "=============================================="
echo " Suit de tests de SynteniToL"
echo "=============================================="

for test_name in "${VALIDATION_TESTS[@]}"; do
    printf '%-40s' "Corriendo ${test_name} ... "
    run_test "$test_name"
    code=$?
    if [[ $code -eq 0 ]]; then
        echo "PASS"
    else
        echo "FAIL (exit ${code})"
        tail -30 "/tmp/syntenitol_test_${test_name}.log"
        FAILED=$((FAILED + 1))
        OVERALL=2
    fi
done

if [[ $SKIP_REGRESSION -eq 1 ]]; then
    echo
    echo "SKIP test_regression.sh (--skip-regression)"
    reg_code=0
else
    echo
    printf '%-40s' "Corriendo test_regression.sh (tarda ~11 min)... "
    set +e
    timeout -k 120 1500 bash tests/test_regression.sh >/tmp/syntenitol_test_regression.log 2>&1
    reg_code=$?
    set -e
    case $reg_code in
        0) echo "PASS" ;;
        1)
            echo "WARNING (resultados distintos a la referencia)"
            tail -20 "/tmp/syntenitol_test_regression.log"
            [[ $OVERALL -eq 0 ]] && OVERALL=1
            ;;
        *)
            echo "FAIL (exit ${reg_code})"
            [[ $reg_code -eq 124 ]] && echo "   (timeout: la regresion supero el limite de 25 min)"
            tail -30 "/tmp/syntenitol_test_regression.log"
            FAILED=$((FAILED + 1))
            OVERALL=2
            ;;
    esac
fi

echo
echo "=============================================="
echo " Resumen"
echo "=============================================="
for entry in "${RESULTS[@]}"; do
    name="${entry%%:*}"
    code="${entry##*:}"
    status="PASS"
    [[ $code -ne 0 ]] && status="FAIL"
    printf '  %-40s %s\n' "${name}" "${status}"
done
if [[ $SKIP_REGRESSION -eq 0 ]]; then
    if [[ $reg_code -eq 0 ]]; then
        printf '  %-40s %s\n' "test_regression.sh" "PASS"
    elif [[ $reg_code -eq 1 ]]; then
        printf '  %-40s %s\n' "test_regression.sh" "WARNING"
    else
        printf '  %-40s %s\n' "test_regression.sh" "FAIL"
    fi
fi
echo

if [[ $OVERALL -eq 0 ]]; then
    echo "Todos los tests pasaron."
elif [[ $OVERALL -eq 1 ]]; then
    echo "Tests pasaron PERO hay diferencias frente a la referencia (revisa el WARNING)."
else
    echo "${FAILED} test(s) fallaron."
fi
exit "$OVERALL"