process merge_mappings {
    publishDir params.outdir, mode: 'copy'
    input:
    path map_files
    output:
    path "all_mappings.map"

    script:
    """
    cat ${map_files} | grep -v "no hits" > all_mappings.map 
    """
}