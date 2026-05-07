process efetch_to_tsv {
    publishDir params.outdir, mode: 'copy'
    label 'retry_backoff'

    input:
    path inputFASTA
    path efetch_script
    output:
    path "protein_to_organism_map.tsv", emit: protein_to_organism_map
    path "assemblies_to_download",      emit: assemblies_to_download
    // path "taxonomy.tsv",                emit: taxonomy_tsv
    path "all_accessions.txt",          emit: all_accessions
    path "no_assembly_list.txt",        emit: no_assembly_list

    script:
    def api_opt = params.ncbi_api_key ? "--ncbi_api_key ${params.ncbi_api_key}" : ""
    """
    echo "Running efetch_to_tsv with inputFASTA: ${inputFASTA}, api_opt: '${api_opt}'"
    conda run -n syntenitol python3 ${efetch_script} \
        --fasta ${inputFASTA} \
        --out_tsv protein_to_organism_map.tsv \
        --out_assemblies assemblies_to_download \
        --out_no_assembly_list no_assembly_list.txt \
        ${api_opt}
    cat assemblies_to_download no_assembly_list.txt > all_accessions.txt
    """
}