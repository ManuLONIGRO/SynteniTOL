process get_taxonomy_info {
    publishDir params.outdir, mode: 'copy'
    label 'retry_backoff'
    input:
        path all_accessions
        path get_taxonomy_info_script
    output:
        path "taxonomy.tsv",   emit: taxonomy_tsv
    
    script:
    def api_opt = params.ncbi_api_key ? "--ncbi_api_key ${params.ncbi_api_key}" : ""
    """
    conda run -n syntenitol python3 ${get_taxonomy_info_script} \
        --input_all_accessions_list ${all_accessions} \
        --output_taxonomy_tsv taxonomy.tsv \
        ${api_opt}
    """

}