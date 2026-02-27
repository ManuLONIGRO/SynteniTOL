#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.help             = params.help ?: false
// params.inputFASTA       = params.inputFASTA ?: null
// params.cogs             = params.cogs ?: null
// params.ncbi_api_key     = params.ncbi_api_key ?: null
//params.color_by_group   = params.color_by_group ?: ""

// params.outdir = params.outdir ?: "results/run_${run_id}"


// Processes
process efetch_to_tsv {
    publishDir params.outdir, mode: 'copy'
    label 'retry_backoff'

    input:
    path inputFASTA
    path efetch_script
    output:
    path "protein_to_organism_map.tsv", emit: protein_to_organism_map
    path "assemblies_to_download",      emit: assemblies_to_download
    path "taxonomy.tsv",                emit: taxonomy_tsv
    path "no_assembly_list.txt",        emit: no_assembly_list

    script:
    def api_opt = params.ncbi_api_key ? "--ncbi_api_key ${params.ncbi_api_key}" : ""
    """
    echo "Running efetch_to_tsv with inputFASTA: ${inputFASTA}, api_opt: '${api_opt}'"
    conda run -n syntenitol python3 ${efetch_script} \
        --fasta ${inputFASTA} \
        --out_tsv protein_to_organism_map.tsv \
        --out_assemblies assemblies_to_download \
        --out_no_assembly_list no_assembly_list.txt \
        --out_taxonomy taxonomy.tsv \
        ${api_opt}
    """
}

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
    if [ -n "${params.ncbi_api_key}" ]; then
        NCBI_API_KEY="${params.ncbi_api_key}"
    fi

    mkdir gbff_assemblies
    printf "%s\n" "${accession_list.join('\n')}" > batch_ids.txt
    rm -f assemblies.zip

    conda run -n syntenitol datasets download genome accession \
        --inputfile batch_ids.txt \
        --include gbff \
        --filename assemblies.zip

    unzip -q assemblies.zip -d gbff_assemblies && rm assemblies.zip    """
}

process download_no_assemblies {
    // Same as download_assemblies but forthe "no_assembly" genomes, which we will download one by one with efetch (since they don't have assemblies, they won't be in the datasets and we have to get them separately)
    maxForks 1
    label 'retry_backoff'
    tag "no_assembly_${task.hash.take(8)}" // tag to identify batches in logs

    input:
        val no_assembly_acc
    output:
        path "gbff_no_assemblies/*.gbff", emit: gbff_no_assembly_file
    script:
    """
    if [ -n "${params.ncbi_api_key}" ]; then
        NCBI_API_KEY="${params.ncbi_api_key}"
    fi

    mkdir -p gbff_no_assemblies
    printf "%s\n" "${no_assembly_acc}" > no_assembly_id.txt

    conda run -n syntenitol efetch -db nuccore -id ${no_assembly_acc} -format gb > gbff_no_assemblies/${no_assembly_acc}.gbff 
    """
}

process rename_gbff_files {
    input:
        tuple path(gbff_dir), path(tsv_file), path(rename_gbff_script)
    output:
        path "genomes"
    
    script:
    """
    conda run -n syntenitol python3 ${rename_gbff_script} \
         --tsv ${tsv_file} \
         --result_zip_ncbi_dir ${gbff_dir}/ncbi_dataset/data \
         --output_dir_genomes genomes
    """
}

process rename_no_assemblies_files {
    input:
    tuple path(gbff_no_assemblies), path(tsv_file), path(rename_no_assemblies_script)
    output:
    path "genomes_no_assemblies"
    script:
    """
    conda run -n syntenitol python3 ${rename_no_assemblies_script} --tsv ${tsv_file} --gbff_no_assemblies ${gbff_no_assemblies} --output_dir_genomes genomes_no_assemblies
    """
}

process merge_genome_directories {
    input:
        val genomes_dirs              // list of absolute paths to each 'genomes' dir
        path genomes_no_assemblies_dir
    output:
        path "genomes_all"

    script:
    def dirs_str = genomes_dirs.collect { it.toString() }.join(' ')
    """
    rm -rf genomes_all
    mkdir -p genomes_all

    for d in ${dirs_str}; do
        cp \$d/*.gbff genomes_all/ || true
    done

    mkdir -p ${genomes_no_assemblies_dir} || true
    cp ${genomes_no_assemblies_dir}/*.gbff genomes_all/ || true
    """
}

process gbff_to_proteomes {
    input:
    tuple path(gbff_files), path(gbff2prot_script)
    output:
    path "proteomes"
    
    script:
    """
    mkdir -p proteomes
    conda run -n syntenitol python3 ${gbff2prot_script} --genomes_dir ${gbff_files} --proteomes_dir proteomes
    """
}

process download_cogs {
    tag { cog }
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

process hmm_build_cogs {
    tag { cog }
    input:
    tuple val(cog), path(cog_file)
    output:
    tuple val(cog), path("profile_${cog}.hmm"), emit: profile_file
    script:
    """
    #mkdir -p cogs_profiles

    base_name="${cog}"
    cdhit -i ${cog_file} -o "\${base_name}_hit" -c 0.9
    awk '/^>/{if(seq)print seq;print;seq="";next}{gsub(/[ \t\r]/,"");seq=seq\$0}END{if(seq)print seq}' "\${base_name}_hit" > "\${base_name}_hit_oneline"
    head -n 1000 "\${base_name}_hit_oneline" > "\${base_name}_hit_oneline_limited"
    mafft --auto --anysymbol "\${base_name}_hit_oneline_limited" > "\${base_name}_alin"
    hmmbuild "profile_\${base_name}.hmm" "\${base_name}_alin"
    #mv "profile_\${base_name}.hmm" "profile_${cog}.hmm"
    """
}

process hmm_search {
    // publishDir "hmm_search_results_dir", mode: 'copy'
    tag { "hmm_search ${cog}" }
    input:
        tuple val(cog), path(profile), val(proteomes_dir)
    output:
    //path "hmm_*.result", emit: hmm_result_file
    tuple val(cog), path("*.result"), emit: hmm_result_file
    
    script:
    """
    profile_name=\$(basename "${profile}" .hmm | sed 's/profile_//')
    
    for proteome in ${proteomes_dir}/*.prot; do
        proteome_name=\$(basename "\$proteome" .prot)
        hmmsearch -E 1e-5 --acc --domtblout "\${proteome_name}_\${profile_name}.result" "${profile}" "\${proteome}"
    done
    """
}

process ids_locustag_mapping {
    input:
    //path hmm_result_file
    //path genomes_dir
    tuple val(cog), path(hmm_result_file), path(genomes_dir), path(ids_locustag_mapping_script)
    output:
    tuple val(cog), path("mappings/*.map"), emit: mapped_file

    script:
    """
    mkdir -p mappings
    base=\$(basename "${hmm_result_file}" .result)
    
    conda run -n syntenitol python3 ${ids_locustag_mapping_script} \
      --genomes_dir ${genomes_dir} \
      --result_file ${hmm_result_file} \
      --output_file mappings/\${base}.map \
      --no_mapped_file mappings/\${base}.no_map
    """
}

process merge_mappings {
    publishDir params.outdir, mode: 'copy'
    input:
    path map_files
    output:
    path "all_mappings.map"

    script:
    """
    cat ${map_files} | grep -v "no hits" > all_mappings.map 
    """
}

process build_dataframe {
    publishDir params.outdir, mode: 'copy'
    input:
    tuple path(mappings_file), path(genomes_dir), path(hmm_results_files), path(map_genes_to_tsv_script)
    output:
    path "dataframe.tsv", emit: dataframe_tsv

    script:
    """
    mkdir -p results_dir
    for f in ${hmm_results_files}; do cp "\$f" results_dir/; done
    conda run -n syntenitol python3 ${map_genes_to_tsv_script} \
      --genomes_dir ${genomes_dir} \
      --mappings_file ${mappings_file} \
      --results_dir results_dir \
      --out_tsv dataframe.tsv
    """
}  

process files_to_synteny {
    publishDir params.outdir, mode: 'copy', pattern: "df_complete.tsv"
    publishDir params.outdir, mode: 'copy', pattern: "candidates.tsv"
	input:
		tuple val(cogs_csv), path(genomes_dir), path(dataframe), path(createFiles2Synteny_script), path(protein_to_organism_map)
	output:
		path "*_presence_binary_data_*", optional: true, emit: presence_binary_data
		path "*_genomic_context_data_*_oriented", optional: true
		path "*_genomic_context_data_*_sector", optional: true, emit: synteny_contexts
        path "df_complete.tsv", emit: df_complete
        path "candidates.tsv", emit: candidates_tsv
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
      --out_candidates_tsv candidates.tsv
	"""
}

process itol_files {
    publishDir params.outdir, mode: 'copy'
    input:
        tuple val(cogs_csv), val(color_by_group), path(genomic_context_file), path(presence_binary_file), path(syntenyTaxOrg_script)
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
      --itol_synteny_file itol_synteny_oriented.txt \
      --synteny_log synteny.log \
      --itol_profiling_file itol_binary.txt
    """
}

process new_format_headers{
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

process taxonomy_itol_files {
    publishDir params.outdir, mode: 'copy', overwrite:true
    input:
    tuple path(taxonomy_tsv_file), path(df_complete_tsv_file), path(fasta_file), path(taxonomy_itol_files_script)
    output:
    path "domain.itol.txt", emit: domain_itol
    path "phylum.itol.txt", emit: phylum_itol
    path "class.itol.txt", emit: class_itol
    script:
    """
    conda run -n syntenitol python3 ${taxonomy_itol_files_script} \
        --taxonomy_tsv ${taxonomy_tsv_file} \
        --df_complete_tsv ${df_complete_tsv_file} \
        --fasta_file ${fasta_file} \
        --domain_itol_file domain.itol.txt \
        --phylum_itol_file phylum.itol.txt \
        --class_itol_file class.itol.txt
    """
}

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

workflow {
    // Create a command.txt with workflow.commandline
    def runDir = file(params.outdir ?: "results/run_${params.run_id}")
    runDir.mkdirs()
    file("${params.outdir}/run_command.txt").text = workflow.commandLine + '\n'

    // Define colors for the messages
    def RED = "\u001B[31m"
    def BLUE = "\u001B[34m"
    def YELLOW = "\u001B[33m"
    def RESET = "\u001B[0m"

    // Print help message if --help is provided
    if (params.help) {
        log.info """
        ${BLUE}SynteniToL Nextflow Pipeline${RESET}
        
        Usage: nextflow run main.nf --inputFASTA <input_fasta_file> --cogs <cog_list> [options]

        ${YELLOW}IMPORTANT:${RESET} the first cog in the list will be used as reference for synteny plots.

        Options:
        --inputFASTA       Input FASTA file with sequences from blastp results (mandatory)
        --cogs             Comma-separated list of COGs (e.g., COG1152,COG1795) (mandatory)
        --color_by_group   Optional parameter to color by group in iTOL files. (e.g. COG1229-COG1029,COG2218,COG2037 COG1229 and COG1029 will be in the same color) 
        --outdir           Output directory (default: results/run_<timestamp>)
        --ncbi_api_key     NCBI API key to increase rate limits (optional but recommended)
        --help             Show this help message and exit
        """
        exit 0
    }

    // Validate mandatory parameters
    if (!params.inputFASTA)   {error "Missing --inputFASTA. Use --help for usage."}
    if (!params.ncbi_api_key) {log.warn "${YELLOW}Warning: No NCBI API key provided. You may encounter rate limits when fetching data from NCBI.${RESET}"}

    // Parse and validate COG list safely (params.cogs may be null/empty)
    def cogs_list = (params.cogs ?: '')
        .tokenize(',')
        .collect { it.trim() }
        .findAll { it }

    if (cogs_list.isEmpty()) { error "Missing --cogs. Use --help for usage." }

    // Check the cogs in the list of cogs. If one COG is other thing that COGXXXX, with XXXX from 0001 to 5950, exit with error. Print the error in red.
    def cog_pattern = ~/^COG(0[0-9]{3}|[1-5][0-9]{3}|5950)$/
    def invalid_cogs = cogs_list.findAll { !(it ==~ cog_pattern) }
    if (invalid_cogs) {
        error "${RED}Invalid COG identifiers found: ${invalid_cogs.join(', ')}. COGs should be in the format COGXXXX, where XXXX is a number from 0001 to 5950.${RESET}"
    }


    // Create channel and call processes with channels
    fasta_ch = channel.fromPath(params.inputFASTA)

    // stage scripts into tasks
    efetch_script_ch = channel.fromPath('bin/efetch_to_tsv.py')
    rename_gbff_script_ch = channel.fromPath('bin/rename_gbff.py')
    gbff2prot_script_ch = channel.fromPath('bin/gbff2prot.py')
    ids_locustag_mapping_script_ch = channel.fromPath('bin/ids_locustag_mapping.py')
    map_genes_to_tsv_script_ch = channel.fromPath('bin/map_genes_to_tsv.py')
    createFiles2Synteny_script_ch = channel.fromPath('bin/createFiles2Synteny.py')
    syntenyTaxOrg_script_ch = channel.fromPath('bin/syntenyTaxOrg.py')
    formatHeaders_script_ch = channel.fromPath('bin/formatHeaders_inputFASTA.py')
    mapping_pid_taxid_script_ch = channel.fromPath('bin/mapping_pid_taxid.py')
    new_formatHeaders_script_ch = channel.fromPath('bin/new_formatHeaders_inputFASTA.py')
    rename_no_assemblies_script_ch = channel.fromPath('bin/rename_no_assemblies.py')
    taxonomy_itol_files_script_ch = channel.fromPath('bin/taxonomy_itol_files.py')
    
    // run efetch_to_tsv process
    //efetch_results = efetch_to_tsv(fasta_ch)
    efetch_results = efetch_to_tsv(fasta_ch,efetch_script_ch)

    // Download assemblies to download in batches of 50 genomes
    efetch_results.assemblies_to_download
        .ifEmpty { error "Assemblies file not found." }
        .splitText()
        .map { it.trim() }
        .filter { it.trim() }
        .toSortedList()
        .flatten()
        .collate(50)
        .set { accession_batches }

    download_assemblies_results = download_assemblies(accession_batches)

    // Download assemblies for those without assembly in NCBI in batches of 50 genomes
    no_assembly_ids = efetch_results.no_assembly_list
                                    .splitText()
                                    .map { it.trim() }
                                    .filter { it }
                                    .toSortedList()
                                    .flatten()
                                    .collate(50)
                                    .set { no_assembly_batches }
    
    no_assemblies_files = download_no_assemblies(no_assembly_batches)
    separate_no_assemblies_files = no_assemblies_files.flatten()

    // Rename no _assemblies_gbff files (this channel may be empty if there are no such genomes)
    rename_no_assemblies_result = rename_no_assemblies_files(
        separate_no_assemblies_files
            .combine(efetch_results.protein_to_organism_map)
            .combine(rename_no_assemblies_script_ch)
    )

    // If there are no "no_assembly" genomes, provide an explicit empty directory placeholder as a *path value* (not a channel)
    rename_no_assemblies_or_empty = rename_no_assemblies_result.ifEmpty { file('empty_no_assemblies') }

    // Rename gbff. download_assemblies_results has just 1 output, so .out doesn't exists.
    genomes_dir = rename_gbff_files(
        download_assemblies_results
            .combine(efetch_results.protein_to_organism_map)
            .combine(rename_gbff_script_ch)
    )

    //genomes_dir.collect().view()
    
    // Merge genomes with (possibly empty) "no_assembly" genomes directory
    all_genomes_dir = merge_genome_directories(
        genomes_dir.collect(),              // → val genomes_dirs (list of work/.../genomes paths)
        rename_no_assemblies_or_empty       // → path genomes_no_assemblies_dir
    )
    //all_genomes_dir.view()
    
    // Creates proteomes from gbff
    proteomes_dir = gbff_to_proteomes(all_genomes_dir.combine(gbff2prot_script_ch))

    // One item per COG for fine-grained caching
    cogs_ch = channel.fromList(cogs_list)

    // Download each COG independently (one task per COG)
    download_cogs_results = download_cogs(cogs_ch)
    cog_files_ch = download_cogs_results.cog_file

    // Build HMM profiles for each COG
    cog_profiles = hmm_build_cogs(cog_files_ch).profile_file
    //cog_profiles.view() //los cogs por separado

    // HMM search of each COG profiles against proteomes
    hmm_search_results = hmm_search(
        cog_profiles
            .combine(proteomes_dir)).hmm_result_file
 
    // Mapping ids_locustag
    // hmm_search emits (cog, [result_files]) per profile; expand to one (cog, result_file) per item
    separate_hmm_results = hmm_search_results.flatMap { cog, files -> files.collect { f -> [cog, f] } }
    
    mapping_hmm_results = ids_locustag_mapping(
        separate_hmm_results
            .combine(all_genomes_dir)
            .combine(ids_locustag_mapping_script_ch)).mapped_file

    final_mapped = merge_mappings(mapping_hmm_results.map { cog, f -> f }.collect())
    
    // Build final dataframe TSV (collect a single flat list of result files)
    map_df = build_dataframe(
        final_mapped
            .combine(all_genomes_dir)
            .combine(separate_hmm_results.map { cog, f -> f }.collect().toList())
            .combine(map_genes_to_tsv_script_ch)
            //.combine(efetch_results.protein_to_organism_map)
    ).dataframe_tsv
    
    // Run synteny file generation on the produced dataframe
    cogs_csv = cogs_list.join(',')

    synteny_results = files_to_synteny(
        channel.value(cogs_csv)
            .combine(all_genomes_dir)
            .combine(map_df)
            .combine(createFiles2Synteny_script_ch)
            .combine(efetch_results.protein_to_organism_map))

    synteny_context_data = synteny_results.synteny_contexts
    presence_binary_data = synteny_results.presence_binary_data
    df_complete_ch = synteny_results.df_complete

    // Build iTOL files
    color_group_opt = params.color_by_group ?: ""
    itol_outputs = itol_files(
        channel.value(cogs_csv)
            .combine(channel.value(color_group_opt))
            .combine(synteny_context_data)
            .combine(presence_binary_data)
            .combine(syntenyTaxOrg_script_ch)
    )

    // format headers of inputFASTA for tree building
    formatted_fasta_ch = new_format_headers(
        fasta_ch
            .combine(all_genomes_dir)
            .combine(new_formatHeaders_script_ch)
    )
    
    // get taxonomy itol files
    taxonomy_itol_files_results = taxonomy_itol_files(
        efetch_results.taxonomy_tsv
            .combine(df_complete_ch)
            .combine(formatted_fasta_ch)
            .combine(taxonomy_itol_files_script_ch)
    )

    // make treefile
    maketree_results = maketree(formatted_fasta_ch)
}