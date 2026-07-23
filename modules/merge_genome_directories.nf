process merge_genome_directories {
    input:
        val genomes_dirs              // list of absolute paths to each 'genomes' dir
        val genomes_no_assemblies_dir
    output:
        path "genomes_all"

    script:
    def dirs_str1 = genomes_dirs.collect { v -> "'${v.toString()}'" }.join(' ')
    def dirs_str2 = genomes_no_assemblies_dir.collect { v -> "'${v.toString()}'" }.join(' ')
    """
    rm -rf genomes_all
    mkdir -p genomes_all

    for d in ${dirs_str1}; do
        cp "\$d"/*.gbff genomes_all/ || true
    done

    for d in ${dirs_str2}; do
        cp "\$d"/*.gbff genomes_all/ || true
    done
    """
}
