process merge_mappings {
    // publishDir params.outdir, mode: 'copy'
    input:
    path map_files, arity: '1..*'
    output:
    path "all_mappings_${task.index}.map"

    script:
    """
    out_file="all_mappings_${task.index}.map"
    : > "\$out_file"

    for f in ${map_files.collect { "'${it}'" }.join(' ')}; do
        # grep returns 1 when no lines match (e.g. file contains only "no hits")
        grep -h -v "no hits" "\$f" >> "\$out_file" || true
    done

    # Remove the repeated lines (a genome processed more than once produced the
    # same mapping line once per copy). The order is preserved.
    awk '!seen[\$0]++' "\$out_file" > "\$out_file".dedup && mv "\$out_file".dedup "\$out_file"
    """
}