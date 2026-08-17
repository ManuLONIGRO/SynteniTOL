process fetch_genbank_records {
    // Same as download_assemblies but forthe "no_assembly" genomes, which we will download one by one with efetch (since they don't have assemblies, they won't be in the datasets and we have to get them separately)
    maxForks 1
    label 'retry_backoff'
    tag "no_assembly_${task.hash.take(8)}" // tag to identify batches in logs

    secret 'NCBI_API_KEY'

    input:
        val no_assembly_acc
    output:
        path "gbff_no_assemblies/*.gbff", emit: gbff_no_assembly_file
    script:
    """
    mkdir -p gbff_no_assemblies

    conda run -n syntenitol efetch \
        -db nuccore \
        -id "${no_assembly_acc}" \
        -format gb \
        > "gbff_no_assemblies/${no_assembly_acc}.gbff"
    """
}