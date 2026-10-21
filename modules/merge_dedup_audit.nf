process merge_dedup_audit {
    publishDir params.outdir, mode: 'copy'
    input:
    val audit_files               // list of per-chunk 'genomes_deduplicated.tsv' files

    output:
    path "genomes_deduplicated.tsv"

    script:
    def audit_list = audit_files.collect { v -> "'${v.toString()}'" }
    def audit_str = audit_list ? audit_list.join(' ') : "''"
    """
    # Merge the per-chunk deduplication audits into a single report. Every
    # chunk file has the same header, so only the first one is kept.
    : > genomes_deduplicated.tsv
    first=1
    for f in ${audit_str}; do
        if [ ! -f "\$f" ]; then continue; fi
        if [ \$first -eq 1 ]; then
            head -n 1 "\$f" >> genomes_deduplicated.tsv
            first=0
        fi
        tail -n +2 "\$f" >> genomes_deduplicated.tsv
    done
    """
}
