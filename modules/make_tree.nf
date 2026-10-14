process make_tree {
    publishDir params.outdir, mode: 'copy'
    input:
    path formatted_fasta
    val base_name
    output:
    path "tree_${base_name}_${params.run_id}.nwk"

    script:
    """
    mafft --anysymbol --auto "${formatted_fasta}" > aligned.fasta
    conda run -n syntenitol bmge -i aligned.fasta -o-of trimmed.fasta -m BLOSUM30 -t AA -h 0.5
    awk '/^>/{n=substr(\$1,2); gsub(/[:(),]/,"_",n); c[n]++; if(c[n]>1) n=n"_"c[n]; print ">"n; next} {print}' trimmed.fasta > unique.fasta
    fasttree unique.fasta > tree_${base_name}_${params.run_id}.nwk
    """
}