# Golden de referencia para el test de regresion

Estos archivos son la salida esperada (golden/referencia) de SynteniToL usando el
FASTA de prueba `data/fwd/small_fwda_test`. El test de regresion re-ejecuta el
comando almacenado en `test_run_command.txt` y compara los archivos `itol_*.txt`
resultantes contra los de este directorio.

## Origen

- Run estable: `results/run_20260916_142328` (16/09/2026).
- Copiados el 17/09/2026 desde ese run.
- Se usa el prefijo `test_` para que estos archivos no caigan bajo el patron
  `itol_*.txt` del `.gitignore` y puedan versionarse en git.
  El comparador mapea `test_itol_X.txt` (golden) -> `itol_X.txt` (run).

## Como correr el test de regresion

    tests/test_regression.sh

o con todos los tests:

    tests/run_all_tests.sh

Si algo cambia, se imprime un WARNING y se escribe un diff detallado en
`<run>/itol_diff_report.txt`. Revisalo para decidir si el cambio corresponde a la
ultima implementacion. Si es un cambio intencional, regenera el golden:

    # 1) corre el pipeline normalmente y localiza el run nuevo
    # 2) reemplaza los archivos de este directorio:
    for f in itol_profiling.txt itol_synteny_conserved.txt \
             itol_synteny_centered_on_fasta_proteins.txt \
             itol_taxonomy_domain.txt itol_taxonomy_phylum.txt itol_taxonomy_class.txt; do
      cp <run_nuevo>/$f test_$f
    done
    cp <run_nuevo>/run_command.txt test_run_command.txt

## Archivos

- `test_run_command.txt` — comando exacto del run estable (se re-ejecuta tal cual).
- `test_itol_*.txt` — los 6 outputs deterministas (sin timestamp) del pipeline.