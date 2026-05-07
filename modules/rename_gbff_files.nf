process rename_gbff_files {
    input:
        tuple path(gbff_dir), path(tsv_file), path(rename_gbff_script)
    output:
        path "genomes"
    
    script:
    """
    conda run -n syntenitol python3 ${rename_gbff_script} \
         --tsv ${tsv_file} \
         --result_zip_ncbi_dir ${gbff_dir}/ncbi_dataset/data \
         --output_dir_genomes genomes
    """
}
