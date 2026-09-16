process generate_synteny_data_fromProtein {
    publishDir params.outdir, mode: 'copy', pattern: "itol_synteny_fasta.txt"

    input:
        tuple val(cogs_csv), path(genomes_dir), path(dataframe), path(createFiles2Synteny_fromProtein_script), val(evalue), path(input_fasta), val(color_by_group)

    output:
        path "itol_synteny_fasta.txt", emit: itol_synteny_fasta

    script:
    """
    out_prefix="synteny_protein"

    extra=""
    if [ -n "${color_by_group}" ]; then
      extra="--color_by_group \"${color_by_group}\""
    fi

    conda run -n syntenitol python3 "${createFiles2Synteny_fromProtein_script}" \
      --input_fasta "${input_fasta}" \
      --input_dir "${genomes_dir}" \
      --dataframe "${dataframe}" \
      --cogs "${cogs_csv}" \
      --evalue ${evalue} \
      \${extra} \
      --itol_synteny_file itol_synteny_fasta.txt
    """
}
