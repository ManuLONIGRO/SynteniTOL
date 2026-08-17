process download_cogs {
    tag { cog }
    maxForks 2
    label 'retry_backoff'

    secret 'NCBI_API_KEY'

    input:
    val cog
    output:
    tuple val(cog), path("cogs/${cog}.fa"),  emit: cog_file
    tuple val(cog), path("cogs/${cog}.tsv"), emit: cog_tsv

    script:
    """
    set -euo pipefail
    mkdir -p cogs

    # Get accession and title for this COG
    efetch -db cdd -id "${cog}" -format docsum \
      | xtract -pattern DocumentSummary -element Accession,Title \
      > "cogs/${cog}.tsv"

    # Download this COG fasta file
    echo "Downloading ${cog} ..."
    wget "https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/fasta/${cog}.fa.gz" -O - \
      | gunzip > "cogs/${cog}.fa"
    """
}
