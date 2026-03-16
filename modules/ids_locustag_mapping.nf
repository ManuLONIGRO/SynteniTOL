process ids_locustag_mapping {
    input:
    //path hmm_result_file
    //path genomes_dir
    tuple val(cog), path(hmm_result_file), path(genomes_dir), path(ids_locustag_mapping_script)
    output:
    tuple val(cog), path("mappings/*.map"), emit: mapped_file

    script:
    """
    mkdir -p mappings
    base=\$(basename "${hmm_result_file}" .result)
    
    conda run -n syntenitol python3 ${ids_locustag_mapping_script} \
      --genomes_dir ${genomes_dir} \
      --result_file ${hmm_result_file} \
      --output_file mappings/\${base}.map \
      --no_mapped_file mappings/\${base}.no_map
    """
}