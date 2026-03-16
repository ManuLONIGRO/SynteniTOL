process itol_files {
    publishDir params.outdir, mode: 'copy'
    input:
        tuple val(cogs_csv), val(color_by_group), path(genomic_context_file), path(presence_binary_file), path(best_goi_tsv) , path(syntenyTaxOrg_script)
    output:
        path "itol_synteny_oriented.txt", emit: itol_synteny
        path "itol_binary.txt", emit: itol_binary
        path "synteny.log", optional: true

    script:
    """
    extra=""
    if [ -n "${color_by_group}" ]; then
      extra="--color_by_group ${color_by_group}"
    fi
    
    conda run -n syntenitol python3 ${syntenyTaxOrg_script} \
      --cogs ${cogs_csv} \
      \${extra} \
      --genomic_context_data ${genomic_context_file} \
      --presence_binary_data ${presence_binary_file} \
      --best_goi_tsv ${best_goi_tsv} \
      --itol_synteny_file itol_synteny_oriented.txt \
      --synteny_log synteny.log \
      --itol_profiling_file itol_binary.txt
    """
}