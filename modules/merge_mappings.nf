process merge_mappings {
    publishDir params.outdir, mode: 'copy'
    input:
    path map_files
    output:
    path "all_mappings.map"

    script:
    """
    : > all_mappings.map

    for f in ${map_files}; do
        grep -h -v "no hits" "\$f" >> all_mappings.map
    done
    """
}