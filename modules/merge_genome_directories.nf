process merge_genome_directories {
    input:
        val genomes_dirs              // list of absolute paths to each 'genomes' dir
        path genomes_no_assemblies_dir
    output:
        path "genomes_all"

    script:
    def dirs_str = genomes_dirs.collect { it.toString() }.join(' ')
    """
    rm -rf genomes_all
    mkdir -p genomes_all

    for d in ${dirs_str}; do
        cp \$d/*.gbff genomes_all/ || true
    done

    mkdir -p ${genomes_no_assemblies_dir} || true
    cp ${genomes_no_assemblies_dir}/*.gbff genomes_all/ || true
    """
}
