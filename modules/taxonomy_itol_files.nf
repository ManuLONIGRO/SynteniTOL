process taxonomy_itol_files {
    publishDir params.outdir, mode: 'copy', overwrite:true
    input:
    tuple path(taxonomy_tsv_file), path(df_complete_tsv_file), path(fasta_file), path(taxonomy_itol_files_script)
    output:
    path "itol_taxonomy_domain.txt", emit: domain_itol
    path "itol_taxonomy_phylum.txt", emit: phylum_itol
    path "itol_taxonomy_class.txt", emit: class_itol
    script:
    """
    conda run -n syntenitol python3 "${taxonomy_itol_files_script}" \
        --taxonomy_tsv "${taxonomy_tsv_file}" \
        --df_complete_tsv "${df_complete_tsv_file}" \
        --fasta_file "${fasta_file}" \
        --domain_itol_file itol_taxonomy_domain.txt \
        --phylum_itol_file itol_taxonomy_phylum.txt \
        --class_itol_file itol_taxonomy_class.txt
    """
}