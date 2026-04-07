process extract_proteomes {
    input:
    tuple path(gbff_files), path(gbff2prot_script)
    output:
    path "proteomes"
    
    script:
    """
    mkdir -p proteomes
    conda run -n syntenitol python3 ${gbff2prot_script} --genomes_dir ${gbff_files} --proteomes_dir proteomes
    """
}