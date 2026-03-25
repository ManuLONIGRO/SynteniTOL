process hmm_search {
    // publishDir "hmm_search_results_dir", mode: 'copy'
    tag "hmm_search ${cog}"
    input:
        tuple val(cog), path(profile), val(proteomes_dir)
    output:
    //path "hmm_*.result", emit: hmm_result_file
    tuple val(cog), path("*.result"), emit: hmm_result_file
    
    script:
    """
    profile_name=\$(basename "${profile}" .hmm | sed 's/profile_//')
    
    for proteome in ${proteomes_dir}/*.prot; do
        proteome_name=\$(basename "\$proteome" .prot)
        hmmsearch -E 1e-5 --acc --domtblout "\${proteome_name}_\${profile_name}.result" "${profile}" "\${proteome}"
    done
    """
}