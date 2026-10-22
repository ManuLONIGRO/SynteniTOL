process get_taxonomy_info {
    publishDir params.outdir, mode: 'copy'
    label 'retry_backoff'

    input:
        path all_accessions
        path get_taxonomy_info_script
    output:
        path "taxonomy.tsv",   emit: taxonomy_tsv
    
    script:
    """
    set -o pipefail
    # edirect prints the full request URL (including api_key=...) to stderr
    # whenever a curl call fails, and Nextflow would store that in
    # .command.log. Merge stderr into stdout and mask the key so it never
    # reaches the task logs, even on transient NCBI errors.
    conda run -n syntenitol python3 "${get_taxonomy_info_script}" \
        --input_all_accessions_list "${all_accessions}" \
        --output_taxonomy_tsv taxonomy.tsv \
        2>&1 | sed -E 's/api_key=[^&" ]*/api_key=***/g'
    """

}