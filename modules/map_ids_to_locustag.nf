process map_ids_to_locustag {
    input:
    //path hmm_result_file
    //path genomes_dir
    tuple val(cog), path(hmm_result_files), path(genomes_dir), path(ids_locustag_mapping_script), val(evalue)
    output:
    tuple val(cog), path("mappings/*.map"), emit: mapped_file

    script:
    """
    mkdir -p mappings
    
    for hmm_result_file in ${hmm_result_files.collect { "'${it}'" }.join(' ')}; do
      base=\$(basename "\${hmm_result_file}" .result)

      conda run -n syntenitol python3 "${ids_locustag_mapping_script}" \
        --genomes_dir "${genomes_dir}" \
        --result_file "\${hmm_result_file}" \
        --gene_name "${cog}" \
        --output_file "mappings/\${base}.map" \
        --no_mapped_file "mappings/\${base}.no_map" \
        --evalue ${evalue}
    done
    """
}