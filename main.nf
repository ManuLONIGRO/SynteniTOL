#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.inputFASTA = file(params.inputFASTA) ?: null  // input file with FASTA sequences
params.cogs = params.cogs ?: null // cogs list

process efetch_to_tsv {
    publishDir "results", mode: 'copy', overwrite: true
    input:
    path inputFASTA
    path efetch_script
    output:
    path "protein_to_organism_map.tsv", emit: protein_to_organism_map
    path "assemblies_to_download", emit: assemblies_to_download
    // path "taxonomy.tsv", emit: taxonomy_tsv
    path "no_assembly_list.txt", emit: no_assembly_list

    script:
    """
    conda run -n syntenitol python3 ${efetch_script} \
        --fasta ${inputFASTA} \
        --out_tsv protein_to_organism_map.tsv \
        --out_assemblies assemblies_to_download \
        --out_no_assembly_list no_assembly_list.txt \
    """
}
// --out_taxonomy taxonomy.tsv \
process mapping_taxid {
    publishDir "results", mode: 'copy', overwrite: true
    input:
    tuple path(taxonomy_tsv), path(protein_to_organism_map), path(mapping_pid_taxid_script)
    output:
    path "protein_taxid.tsv", emit: protein_taxid_map
    script:
    """
    conda run -n syntenitol python3 ${mapping_pid_taxid_script} \
        --protein_to_organism_map ${protein_to_organism_map} \
        --taxonomy_tsv ${taxonomy_tsv} \
        --out_pid_taxid_tsv protein_taxid.tsv
    """
}
process download_assemblies {
    // errorStrategy 'retry'
    // maxRetries 1   // second attempt if the download fails (e.g. network hiccup)
    input:
    path assemblies_to_download
    output:
    path "gbff_assemblies"
    
    script:
    """
    mkdir gbff_assemblies
    rm -f assemblies.zip
    conda run -n syntenitol datasets download genome accession --inputfile ${assemblies_to_download} --include gbff --filename assemblies.zip
    unzip assemblies.zip -d gbff_assemblies
    """
}
process download_no_assemblies {
    input:
    val no_assembly_acc
    output:
    path "gbff_no_assemblies/*.gbff", emit: gbff_no_assembly_file
    script:
    """ 
    #set -euo pipefail
    mkdir -p gbff_no_assemblies
    conda run -n syntenitol efetch -db nuccore -id ${no_assembly_acc} -format gb > gbff_no_assemblies/${no_assembly_acc}.gbff 
    """
}
process rename_gbff_files {
    input:
    path gbff_dir 
    path tsv_file
    path rename_gbff_script
    output:
    path "genomes"
    
    script:
    """
    conda run -n syntenitol python3 ${rename_gbff_script} --tsv ${tsv_file} --result_zip_ncbi_dir ${gbff_dir}/ncbi_dataset/data --output_dir_genomes genomes
    """
}

process rename_no_assemblies_files {
    input:
    path gbff_no_assemblies
    path tsv_file
    path rename_no_assemblies_script
    output:
    path "genomes_no_assemblies"
    script:
    """
    conda run -n syntenitol python3 ${rename_no_assemblies_script} --tsv ${tsv_file} --gbff_no_assemblies ${gbff_no_assemblies} --output_dir_genomes genomes_no_assemblies
    """
}

process merge_genome_directories {
    input:
        tuple path(genomes_dir), path(genomes_no_assemblies_dir)
    output:
        path "genomes_all"

    script:
    """
    rm -rf genomes_all
    mkdir -p genomes_all
    cp ${genomes_dir}/*.gbff genomes_all/ || true
    mkdir -p ${genomes_no_assemblies_dir} || true
    cp ${genomes_no_assemblies_dir}/*.gbff genomes_all/ || true
    """
}

process gbff_to_proteomes {
    input:
    path gbff_files
    path gbff2prot_script
    output:
    path "proteomes"
    
    script:
    """
    mkdir -p proteomes
    conda run -n syntenitol python3 ${gbff2prot_script} --genomes_dir ${gbff_files} --proteomes_dir proteomes
    """
}
//conda run -n syntenitol python3 ${gbff2prot_script} --genomes_dir ${gbff_no_assemblies_files} --proteomes_dir proteomes
//#conda run -n syntenitol python3 gbff2prot.py --genomes_dir ${gbff_files} --proteomes_dir proteomes
process download_cogs {
    input:
    val cogs_list // List of cogs from the user. COG1152,COG1795,...
    output:
    path "cogs/*.fa", emit: cog_file
    path "cogs.tsv", emit: cogs_tsv
    
    script:
    """
    set -euo pipefail
    mkdir -p cogs
    
    # Creates a file with the list of cogs
    printf "%s\n" ${cogs_list.collect{ cog -> "'${cog}'" }.join(' ')} > cogs_ids.txt

    # Get accesion and title for each COG
    efetch -db cdd -id ${cogs_list.join(',')} -format docsum | xtract -pattern DocumentSummary -element Accession,Title > cogs.tsv
    
    # Download each COG fasta file
    while read cog; do
        echo "Downloading \$cog ..."
        wget "https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/fasta/\$cog.fa.gz" -O - | gunzip > cogs/\$cog.fa
    done < cogs_ids.txt
    """
}

process hmm_build_cogs {
    tag "hmm_build on ${cog_file}"
    input:
    path cog_file
    output:
    path "profile_*.hmm",emit:profile_file
    script:
    """
    #mkdir -p cogs_profiles

    base_name=\$(basename "${cog_file}" .fa)
    cdhit -i ${cog_file} -o "\${base_name}_hit" -c 0.9
    mafft --auto --anysymbol "\${base_name}_hit" > "\${base_name}_alin"
    hmmbuild "profile_\${base_name}.hmm" "\${base_name}_alin"
    #mv "profile_\${base_name}.hmm" cogs_profiles/
    """
}

process hmm_search {
    tag "hmm_search ${profile.simpleName}"
    input:
        tuple path(profile), val(proteomes_dir)
    output:
    //path "hmm_*.result", emit: hmm_result_file
    path "*.result", emit: hmm_result_file
    
    script:
    """
    profile_name=\$(basename "${profile}" .hmm | sed 's/profile_//')
    
    for proteome in ${proteomes_dir}/*.prot; do
        proteome_name=\$(basename "\$proteome" .prot)
        hmmsearch -E 1e-5 --acc --domtblout "\${proteome_name}_\${profile_name}.result" "${profile}" "\${proteome}"
    done
    """
}
// #mkdir -p hmm
//#for file in ${profile}/*.hmm; do
    //#hmmsearch "${profile}" "\${proteome}" > "hmm_\${proteome_name}_\${profile_name}.result"
    //#mv "\${proteome_name}_\$profile_name.result" hmm_search_results_dir/

process ids_locustag_mapping {
    input:
    //path hmm_result_file
    //path genomes_dir
    tuple path(hmm_result_file), val(genomes_dir), path(ids_locustag_mapping_script)
    output:
    path "mappings/*.map", emit: mapped_file

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
/*
#--metrics_file mappings/\${base}.metrics.tsv
    # ensure files exist even if empty
    # touch mappings/\${base}.mapped.txt mappings/\${base}.no_mapped.txt mappings/\${base}.metrics.tsv
*/
//mapped lines (protein_id -> locus_tag)
    //conda run -n syntenitol python3 ids_locustag_mapping.py \
process merge_mappings {
    publishDir "results", mode: 'copy', overwrite: true
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
    publishDir "results", mode: 'copy', overwrite: true
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
//  conda run -n syntenitol python3 map_genes_to_tsv.py \  

    //, path(protein_map)    

process files_to_synteny {
    publishDir "results", mode: 'copy', pattern: "df_complete.tsv", overwrite: true
	input:
		tuple val(cogs_csv), path(genomes_dir), path(dataframe), path(createFiles2Synteny_script), path(protein_to_organism_map)
	output:
		path "*_presence_binary_data_*", optional: true, emit: presence_binary_data
		path "*_genomic_context_data_*_oriented", optional: true
		path "*_genomic_context_data_*_sector", optional: true, emit: synteny_contexts
        path "df_complete.tsv", emit: df_complete
	script:
	"""
	out_prefix="synteny"
	
    conda run -n syntenitol python3 ${createFiles2Synteny_script} \
	  --input_dir ${genomes_dir} \
	  --cogs ${cogs_csv} \
	  --dataframe ${dataframe} \
	  --out_prefix \${out_prefix} \
      --protein_to_organism_map_tsv ${protein_to_organism_map} \
      --out_tsv df_complete.tsv
	"""
}
//#conda run -n syntenitol python3 createFiles2Synteny.py \
process itol_files {
    publishDir "results", mode: 'copy', overwrite: true
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
    publishDir "results", mode: 'copy', overwrite:true
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

// process format_headers_fasta {
//     publishDir "results", mode: 'copy', overwrite: true
//     input:
//     tuple path(inputFASTA), path(df_complete), path(formatHeaders_script_ch)

//     output:
//     path "formatted_headers.fasta"

//     script:
//     """
//     conda run -n syntenitol python3 ${formatHeaders_script_ch} \
//         --input_fasta ${inputFASTA} \
//         --dataframe ${df_complete} \
//         --output_fasta formatted_headers.fasta
//     """
// }
//  path(protein_to_organism_map),
// --protein_to_organism_map ${protein_to_organism_map} \
//#conda run -n syntenitol python3 syntenyTaxOrg.py \
process maketree {
    publishDir "results", mode: 'copy', overwrite: true
    input:
    path formatted_fasta
    output:
    path "tree.nwk"

    script:
    """
    mafft --auto ${formatted_fasta} > aligned.fasta
    conda run -n syntenitol bmge -i aligned.fasta -o-of trimmed.fasta -m BLOSUM30 -t AA -h 0.5
    fasttree trimmed.fasta > tree.nwk
    """
}

workflow {
    
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
    // run efetch_to_tsv process
    //efetch_results = efetch_to_tsv(fasta_ch)
    efetch_results = efetch_to_tsv(fasta_ch,efetch_script_ch)
    
    // Creating pid_taxid tsv No tiene sentido porque tiene los taxid del lineage pero no del organismo.""
    // protein_taxid_map_tsv = mapping_taxid(efetch_results.taxonomy_tsv
    //     .combine(efetch_results.protein_to_organism_map)
    //     .combine(mapping_pid_taxid_script_ch))

    // Download assemblies
    download_assemblies_results = download_assemblies(efetch_results.assemblies_to_download)
    // Download assemblies for those without assembly in NCBI
    no_assembly_ids = efetch_results.no_assembly_list
                                    .splitText()
                                    .map { it.trim() }
                                    .filter { it }
    no_assemblies_files = download_no_assemblies(no_assembly_ids)
    separate_no_assemblies_files = no_assemblies_files.flatten()
    // Rename no _assemblies_gbff files (this channel may be empty if there are no such genomes)
    rename_no_assemblies_result = rename_no_assemblies_files(separate_no_assemblies_files, efetch_results.protein_to_organism_map, rename_no_assemblies_script_ch)
    // If there are no "no_assembly" genomes, provide an explicit empty directory placeholder as a *path value* (not a channel)
    rename_no_assemblies_or_empty = rename_no_assemblies_result.ifEmpty { file('empty_no_assemblies') }

    // Rename gbff. download_assemblies_results has just 1 output, so .out doesn't exists.
    genomes_dir = rename_gbff_files(download_assemblies_results, efetch_results.protein_to_organism_map,rename_gbff_script_ch)
    // Merge genomes with (possibly empty) "no_assembly" genomes directory
    all_genomes_dir = merge_genome_directories(genomes_dir.combine(rename_no_assemblies_or_empty))

    // Creates proteomes from gbff
    proteomes_dir = gbff_to_proteomes(all_genomes_dir,gbff2prot_script_ch)

    // convert string "COG1152,COG1795" to list ['COG1152','COG1795']
    cogs_list = params.cogs.tokenize(',')*.trim()
    //cog_list_ch = channel.fromList(cogs_list) //COG1121\nCOG1795
    //cog_list_ch.view()    
    // Download COGs from the list in the input cog
    download_cogs_results = download_cogs(cogs_list)
    cog_files_ch = download_cogs_results.cog_file

    // Flatten channel of cog files
    separate_cogs_ch = cog_files_ch.flatten()

    // Build HMM profiles for each COG
    cog_profiles = hmm_build_cogs(separate_cogs_ch).profile_file
    //cog_profiles.view() //los cogs por separado

    // HMM search of each COG profiles against proteomes
    hmm_search_results = hmm_search(cog_profiles.combine(proteomes_dir)).hmm_result_file
 
    // Mapping ids_locustag
    separate_hmm_results = hmm_search_results.flatten()
    //separate_hmm_results.view()
    mapping_hmm_results = ids_locustag_mapping(separate_hmm_results.combine(all_genomes_dir).combine(ids_locustag_mapping_script_ch)).mapped_file
    final_mapped = merge_mappings(mapping_hmm_results.collect())
    
    // Build final dataframe TSV (collect a single flat list of result files)
    map_df = build_dataframe(
        final_mapped
            .combine(all_genomes_dir)
            .combine(separate_hmm_results.collect().toList())
            .combine(map_genes_to_tsv_script_ch)
            //.combine(efetch_results.protein_to_organism_map)
    ).dataframe_tsv
    
    // Run synteny file generation on the produced dataframe
    cogs_csv = params.cogs
    synteny_results = files_to_synteny(channel.value(cogs_csv)
                                .combine(all_genomes_dir)
                                .combine(map_df)
                                .combine(createFiles2Synteny_script_ch)
                                .combine(efetch_results.protein_to_organism_map))
    synteny_context_data = synteny_results.synteny_contexts
    presence_binary_data = synteny_results.presence_binary_data
    df_complete = synteny_results.df_complete

    // Build iTOL files
    color_group_opt = params.color_by_group ?: ""
    itol_outputs = itol_files(channel.value(cogs_csv)
                                .combine(channel.value(color_group_opt))
                                .combine(synteny_context_data)
                                .combine(presence_binary_data)
                                .combine(syntenyTaxOrg_script_ch))
    // format headers of inputFASTA for tree building
    formatted_fasta_ch = new_format_headers(fasta_ch
                                                .combine(all_genomes_dir)
                                                .combine(new_formatHeaders_script_ch))
    
    // make treefile
    maketree_results = maketree(formatted_fasta_ch)
}


// Comento los dos procesos para generar phylum.itol.txt
// El input deberia ser el resultado de eliminar a los organismos que solo tengan genero y no especie (awk)
/*
process extract_headers {
    tag "$fasta"
    input:
        path fasta
    output:
        path "headers.txt"

    script:
    """
    grep ">" $fasta | sed 's/>//' > headers.txt
    """
}

process process_headers {
    // publish results
    publishDir "results", mode: 'copy'
    
    input:
        path headers
        path itol_script
    output:
        path "itol_results"

    script:
    """
    python3 $itol_script $headers itol_results/
    """
}

// Mandatory validations


if (!params.cogs) {
    log.error "❌ Mandatory parameter. Please provide a COGs list with --cogs, like --cogs COG1152,COG1795"
}
if (!params.inputFASTA) {
    log.error "❌ Mandatory parameter. Please provide the result of blastp with --inputFASTA, like --inputFASTA result_from_blastp.fasta"
}
*/