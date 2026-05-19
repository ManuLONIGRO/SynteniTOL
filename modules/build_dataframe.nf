process build_dataframe {
    publishDir params.outdir, mode: 'copy'
    input:
    path mappings_file
    path genomes_dir
    path hmm_results_files, arity: '1..*'
    path map_genes_to_tsv_script
    output:
    path "dataframe.tsv", emit: dataframe_tsv

    script:
    """
    mkdir -p results_dir
    for f in ${hmm_results_files}; do cp "\$f" results_dir/; done
    conda run -n syntenitol python3 ${map_genes_to_tsv_script} \
      --genomes_dir ${genomes_dir} \
      --mappings_file ${mappings_file} \
      --results_dir results_dir \
      --out_tsv dataframe.tsv
    """
}