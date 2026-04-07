process reformat_fasta_headers {
    publishDir params.outdir, mode: 'copy', overwrite:true
    input:
    tuple path(inputFASTA), path(genomes_dir), path(new_formatHeaders_script)
    output:
    path "formatted_headers.fasta"
    script:
    """
    conda run -n syntenitol python3 ${new_formatHeaders_script} \
        --input_fasta ${inputFASTA} \
        --genomes_dir ${genomes_dir} \
        --output_fasta formatted_headers.fasta
    """
}