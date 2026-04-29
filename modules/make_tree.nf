process make_tree {
    publishDir params.outdir, mode: 'copy'
    input:
    path formatted_fasta
    
    script:
    // Extract filename without extension from inputFASTA parameter
    def input_file = file(params.inputFASTA)
    def base_name = input_file.baseName
    def tree_name = "${base_name}_${params.run_id}"
    
    output:
    path "tree_${tree_name}.nwk"
    
    """
    mafft --auto ${formatted_fasta} > aligned.fasta
    conda run -n syntenitol bmge -i aligned.fasta -o-of trimmed.fasta -m BLOSUM30 -t AA -h 0.5
    fasttree trimmed.fasta > tree_${tree_name}.nwk
    """
}