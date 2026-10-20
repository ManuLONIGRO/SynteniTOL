process download_assemblies {
    maxForks 1
    label 'retry_backoff'
    tag "batch_${task.hash.take(8)}" // tag to identify batches in logs

    input:
    val accession_list
    output:
    path "gbff_assemblies"
    
    script:
    
    """
    mkdir gbff_assemblies
    printf "%s\n" "${accession_list.join('\n')}" > batch_ids.txt
    rm -f assemblies.zip

    conda run -n syntenitol datasets download genome accession \
        --inputfile batch_ids.txt \
        --include gbff \
        --filename assemblies.zip

    unzip -q assemblies.zip -d gbff_assemblies && rm assemblies.zip    """
}
