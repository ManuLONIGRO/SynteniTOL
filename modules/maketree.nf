process maketree {
    publishDir params.outdir, mode: 'copy'
    input:
    path formatted_fasta
    output:
    path "tree_${params.run_id}.nwk"

    script:
    """
    mafft --auto ${formatted_fasta} > aligned.fasta
    conda run -n syntenitol bmge -i aligned.fasta -o-of trimmed.fasta -m BLOSUM30 -t AA -h 0.5
    fasttree trimmed.fasta > tree_${params.run_id}.nwk
    """
}