process generate_synteny_data {
    publishDir params.outdir, mode: 'copy', pattern: "df_complete.tsv"
    publishDir params.outdir, mode: 'copy', pattern: "candidates.tsv"
	publishDir params.outdir, mode: 'copy', pattern: "best_goi.tsv"
    input:
		tuple val(cogs_csv), path(genomes_dir), path(dataframe), path(createFiles2Synteny_script), path(protein_to_organism_map)
	output:
		path "*_presence_binary_data_*", optional: true, emit: presence_binary_data
		path "*_genomic_context_data_*_oriented", optional: true
		path "*_genomic_context_data_*_sector", optional: true, emit: synteny_contexts
        path "df_complete.tsv", emit: df_complete
        path "candidates.tsv", emit: candidates_tsv
        path "best_goi.tsv", emit: best_goi_tsv

	script:
	"""
	out_prefix="synteny"
	
    conda run -n syntenitol python3 ${createFiles2Synteny_script} \
	  --input_dir ${genomes_dir} \
	  --cogs ${cogs_csv} \
	  --dataframe ${dataframe} \
	  --out_prefix \${out_prefix} \
      --protein_to_organism_map_tsv ${protein_to_organism_map} \
      --out_tsv df_complete.tsv \
      --out_candidates_tsv candidates.tsv \
      --out_best_goi_tsv best_goi.tsv
	"""
}