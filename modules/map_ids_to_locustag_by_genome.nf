process map_ids_to_locustag_by_genome {
    tag "map_ids_to_locustag ${genome}"

    input:
        tuple val(genome), val(cogs), path(result_files), path(genomes_dir), path(script), val(evalue)

    output:
        tuple val(genome), path("mappings/*.map"), emit: mapped_file, optional: true

    script:
    // Normalise to lists (defensive: a group should always be a list).
    def cogList  = cogs instanceof List ? cogs : [cogs]
    def fileList = result_files instanceof List ? result_files : [result_files]
    def files_arg = fileList.collect { f -> "'${f}'" }.join(' ')
    def cogs_arg  = cogList.join(',')
    """
    mkdir -p mappings

    conda run -n syntenitol python3 "${script}" \
        --genomes_dir "${genomes_dir}" \
        --genome "${genome}" \
        --result_files ${files_arg} \
        --gene_names "${cogs_arg}" \
        --out_dir mappings \
        --evalue ${evalue}
    """
}
