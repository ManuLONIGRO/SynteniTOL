process rename_no_assemblies_files {
    input:
    tuple path(gbff_no_assemblies), path(tsv_file), path(rename_no_assemblies_script)
    output:
    path "genomes_no_assemblies"
    script:
    """
    conda run -n syntenitol python3 ${rename_no_assemblies_script} --tsv ${tsv_file} --gbff_no_assemblies ${gbff_no_assemblies} --output_dir_genomes genomes_no_assemblies
    """
}