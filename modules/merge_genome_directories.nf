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

    # The genomes are HARD LINKED instead of copied: the bytes of every
    # assembly are stored only once in work/ (the downloaded assembly, the
    # renamed genome and genomes_all all point to the same inode), so the
    # merged directory costs almost no extra disk space.
    # If the hard link cannot be created (e.g. different filesystems) a
    # normal copy is used as a fallback, so the pipeline always works.
    n_linked=0
    n_copied=0

    merge_dir() {
        for f in "\$1"/*.gbff; do
            [ -e "\$f" ] || continue
            base=\$(basename "\$f")
            # readlink -f: link the real file, not a symlink pointing to it.
            src=\$(readlink -f "\$f")
            [ -n "\$src" ] || continue
            if ln -f "\$src" "genomes_all/\$base" 2>/dev/null; then
                n_linked=\$((n_linked + 1))
            else
                # Fallback for filesystems without hard link support: the copy
                # goes through a temporary name and is renamed, so an existing
                # destination is never truncated in place (it could share its
                # inode with another genome).
                tmp="genomes_all/.cp_tmp.\$\$"
                if cp "\$src" "\$tmp" 2>/dev/null && mv -f "\$tmp" "genomes_all/\$base"; then
                    n_copied=\$((n_copied + 1))
                else
                    rm -f "\$tmp"
                fi
            fi
        done
    }

    if [ -n "${dirs_str1}" ]; then
        for d in ${dirs_str1}; do
            merge_dir "\$d"
        done
    fi

    if [ -n "${dirs_str2}" ]; then
        for d in ${dirs_str2}; do
            merge_dir "\$d"
        done
    fi

    echo "[=] genomes_all: \${n_linked} hard links | \${n_copied} copies"
    """
}
