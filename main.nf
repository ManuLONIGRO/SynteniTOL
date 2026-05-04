#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES - PROCESSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { efetch_to_tsv                 } from './modules/efetch_to_tsv.nf'
include { get_taxonomy_info             } from './modules/get_taxonomy_info.nf'
include { download_assemblies           } from './modules/download_assemblies.nf'
include { fetch_genbank_records         } from './modules/fetch_genbank_records.nf'
include { rename_gbff_files             } from './modules/rename_gbff_files.nf'
include { rename_genbank_files          } from './modules/rename_genbank_files.nf'
include { merge_genome_directories      } from './modules/merge_genome_directories.nf'
include { extract_proteomes             } from './modules/extract_proteomes.nf'
include { download_cogs                 } from './modules/download_cogs.nf'
include { hmm_build_cogs                } from './modules/hmm_build_cogs.nf'
include { hmm_search                    } from './modules/hmm_search.nf'
include { map_ids_to_locustag           } from './modules/map_ids_to_locustag.nf'
include { merge_mappings                } from './modules/merge_mappings.nf'
include { build_dataframe               } from './modules/build_dataframe.nf'
include { generate_synteny_data         } from './modules/generate_synteny_data.nf'
include { generate_itol_files           } from './modules/generate_itol_files.nf'
include { reformat_fasta_headers        } from './modules/reformat_fasta_headers.nf'
include { taxonomy_itol_files           } from './modules/taxonomy_itol_files.nf'
include { make_tree                     } from './modules/make_tree.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Create a command.txt with workflow.commandline in the output directory for reproducibility
    def runDir = file(params.outdir ?: "results/run_${params.run_id}")
    runDir.mkdirs()
    file("${params.outdir}/run_command.txt").text = workflow.commandLine + '\n'

    // Define colors for the messages
    def RED     = "\u001B[31m"
    def BLUE    = "\u001B[34m"
    def YELLOW  = "\u001B[33m"
    def RESET   = "\u001B[0m"

    // Print help message if --help is provided
    if (params.help) {
        log.info """
        ${BLUE}SynteniToL Nextflow Pipeline${RESET}
        
        Usage: nextflow run main.nf --inputFASTA <input_fasta_file> --cogs <cog_list> [options]

        ${YELLOW}IMPORTANT:${RESET} the first profile in the combined list (COGs + custom profiles) is used as reference for synteny plots.

        Options:
        --inputFASTA       Input FASTA file with sequences from blastp results (mandatory)
        --cogs             Comma-separated list of COGs (e.g., COG1152,COG1795) (mandatory)
        --custom_hmm_profiles Optional comma-separated list of custom HMM profile paths to include (e.g., profile_x.hmm,profile_y.hmm)
        --query_profile    Optional profile label to use as reference for inputFASTA (must exist in COGs/custom profiles)
        --color_by_group   Optional parameter to color by group in iTOL files. (e.g. COG1229-COG1029,COG2218,COG2037 COG1229 and COG1029 will be in the same color) 
        --evalue           E-value threshold used in search/filtering steps (default: 1e-5)
        --outdir           Output directory (default: results/run_<timestamp>)
        --ncbi_api_key     NCBI API key to increase rate limits (optional but recommended)
        --help             Show this help message and exit
        """
        exit 0
    }


    /////////////////////////////////////////////////////////////
    /*                 Validate parameters                     */
    /////////////////////////////////////////////////////////////

    // Validate mandatory parameters
    if (!params.inputFASTA)   {error "Missing --inputFASTA. Use --help for usage."}
    def inputFile = file(params.inputFASTA)
    if (!inputFile.exists()) {
        error "The input file ${params.inputFASTA} doesn't exist."
    }

    if (!params.ncbi_api_key) {log.warn "${YELLOW}Warning: No NCBI API key provided. You may encounter rate limits when fetching data from NCBI.${RESET}"}

    // Get the name of the input file without the extension for the tree name file.
    def base_name = java.nio.file.Paths.get(params.inputFASTA.toString()).fileName.toString().replaceFirst(/\.[^.]+$/, '')

    // Define evalue parameter channel
    def evalue = channel.value(params.evalue)
    // Parse and validate COG list safely (params.cogs may be null/empty)
    def cogs_list = (params.cogs ?: '')
        .tokenize(',')
        .collect { v -> v.trim() }
        .findAll { v -> v }
    if (cogs_list.isEmpty()) { error "Missing --cogs. Use --help for usage." }

    // Optional custom HMM profiles provided by user (comma-separated paths)
    def custom_profiles_list = (params.custom_hmm_profiles ?: '')
        .tokenize(',')
        .collect { v -> v.trim() }
        .findAll { v -> v }
    def custom_profile_tuples = custom_profiles_list.collect { profilePath ->
        def profileFile = file(profilePath)
        if (!profileFile.exists()) {
            error "Custom HMM profile not found: ${profilePath}"
        }
        def profileName = java.nio.file.Paths.get(profilePath).fileName.toString().replaceFirst(/\.hmm$/, '').replaceFirst(/^profile_/, '')
        [profileName, profileFile]
    }
    def custom_profile_names = custom_profile_tuples.collect { entry -> entry[0] }

    // Check the cogs in the list of cogs. If one COG is other thing that COGXXXX, with XXXX from 0001 to 5950, exit with error. Print the error in red.
    def cog_pattern = ~/^COG(0[0-9]{3}|[1-5][0-9]{3}|5950)$/
    def invalid_cogs = cogs_list.findAll { cog -> !(cog ==~ cog_pattern) }
    if (invalid_cogs) {
        error "${RED}Invalid COG identifiers found: ${invalid_cogs.join(', ')}. COGs should be in the format COGXXXX, where XXXX is a number from 0001 to 5950.${RESET}"
    }


    /////////////////////////////////////////////////////////////
    /*               Channels for processes                    */
    /////////////////////////////////////////////////////////////    

    // Create channel and call processes with channels
    fasta_ch                            = channel.fromPath(params.inputFASTA)
    efetch_script_ch                    = channel.fromPath('bin/efetch_to_tsv.py')
    get_taxonomy_info_script_ch         = channel.fromPath('bin/get_taxonomy_info.py')
    rename_gbff_script_ch               = channel.fromPath('bin/rename_gbff.py')
    gbff2prot_script_ch                 = channel.fromPath('bin/gbff2prot.py')
    ids_locustag_mapping_script_ch      = channel.fromPath('bin/ids_locustag_mapping.py')
    map_genes_to_tsv_script_ch          = channel.fromPath('bin/map_genes_to_tsv.py')
    createFiles2Synteny_script_ch       = channel.fromPath('bin/createFiles2Synteny.py')
    syntenyTaxOrg_script_ch             = channel.fromPath('bin/syntenyTaxOrg.py')
    new_formatHeaders_script_ch         = channel.fromPath('bin/new_formatHeaders_inputFASTA.py')
    rename_no_assemblies_script_ch      = channel.fromPath('bin/rename_no_assemblies.py')
    taxonomy_itol_files_script_ch       = channel.fromPath('bin/taxonomy_itol_files.py')


    /////////////////////////////////////////////////////////////
    /*                      Run processes                      */
    ///////////////////////////////////////////////////////////// 

    efetch_results                      = efetch_to_tsv(
                                            fasta_ch,
                                            efetch_script_ch
                                            )
    
    
    taxonomy_info_results               = get_taxonomy_info(
                                            efetch_results.all_accessions,
                                            get_taxonomy_info_script_ch
                                        )

    // Download assemblies to download in batches of 50 genomes to not overload NCBI servers
    efetch_results.assemblies_to_download
        .ifEmpty { error "Assemblies file not found." }
        .splitText()
        .map { v -> v.trim() }
        .filter { v -> v.trim() }
        .toSortedList()
        .flatten()
        .collate(50)
        .set { accession_batches }

    download_assemblies_results     = download_assemblies(accession_batches)

    /*
    no_assembly_ids                 = efetch_results.no_assembly_list
                                        .splitText()
                                        .map { v -> v.trim() }
                                        .filter { v -> v }

    no_assemblies_files             = fetch_genbank_records(no_assembly_ids)
    separate_no_assemblies_files    = no_assemblies_files.flatten()

    

    rename_no_assemblies_result     = rename_genbank_files(
                                        separate_no_assemblies_files
                                            .combine(efetch_results.protein_to_organism_map)
                                            .combine(rename_no_assemblies_script_ch)
                                    )
    
    // If there are no "no_assembly" genomes, provide an explicit empty directory placeholder as a *path value* (not a channel)
    rename_no_assemblies_or_empty   = rename_no_assemblies_result.ifEmpty { file('empty_no_assemblies') }
    */
    rename_no_assemblies_or_empty = file('empty_no_assemblies')
    
    genomes_dir                 = rename_gbff_files(
                                    download_assemblies_results
                                        .combine(efetch_results.protein_to_organism_map)
                                        .combine(rename_gbff_script_ch)
                                )
    
    

    all_genomes_dir             = merge_genome_directories(
                                    genomes_dir.collect(),              // val genomes_dirs (list of work/.../genomes paths)
                                    rename_no_assemblies_or_empty.collect()       // val genomes_no_assemblies_dir
                                )

    proteomes_dir               = extract_proteomes(all_genomes_dir.combine(gbff2prot_script_ch))

    
    cogs_ch                     = channel.fromList(cogs_list)
    download_cogs_results       = download_cogs(cogs_ch)
    cog_files_ch                = download_cogs_results.cog_file
    cog_profiles                = hmm_build_cogs(cog_files_ch).profile_file
    custom_profiles_ch          = custom_profile_tuples ? channel.fromList(custom_profile_tuples) : channel.empty()
    all_profiles_ch             = cog_profiles.mix(custom_profiles_ch)
    hmm_search_results          = hmm_search(
                                    all_profiles_ch
                                        .combine(proteomes_dir)
                                        .combine(evalue)).hmm_result_file
    separate_hmm_results        = hmm_search_results.flatMap { cog, files -> files.collect { f -> [cog, f] } }

    
    mapping_hmm_results         = map_ids_to_locustag(
                                    hmm_search_results
                                        .combine(all_genomes_dir)
                                        .combine(ids_locustag_mapping_script_ch)
                                        .combine(evalue)).mapped_file
    final_mapped                = merge_mappings(mapping_hmm_results.map { _cog, f -> f }.collect())


    map_df                      = build_dataframe(
                                    final_mapped
                                        .combine(all_genomes_dir)
                                        .combine(separate_hmm_results.map { _cog, f -> f }.collect().toList())
                                        .combine(map_genes_to_tsv_script_ch)
                                ).dataframe_tsv


    all_gene_labels             = cogs_list + custom_profile_names
    def query_profile = params.query_profile?.toString()?.trim()
    if (query_profile) {
        if (!(query_profile in all_gene_labels)) {
            error "Invalid --query_profile '${query_profile}'. It must match one of: ${all_gene_labels.join(', ')}"
        }
    } else if (custom_profile_names) {
        query_profile = custom_profile_names[0]
    } else {
        query_profile = cogs_list[0]
    }
    ordered_gene_labels         = [query_profile] + all_gene_labels.findAll { gene -> gene != query_profile }
    cogs_csv                    = ordered_gene_labels.join(',')
    synteny_results             = generate_synteny_data(
                                    channel.value(cogs_csv)
                                        .combine(all_genomes_dir)
                                        .combine(map_df)
                                        .combine(createFiles2Synteny_script_ch)
                                        .combine(efetch_results.protein_to_organism_map)
                                        .combine(evalue))


    synteny_context_data        = synteny_results.synteny_contexts
    presence_binary_data        = synteny_results.presence_binary_data
    df_complete_ch              = synteny_results.df_complete
    best_goi_tsv_ch             = synteny_results.best_goi_tsv


    // Build iTOL files
    color_group_opt             = params.color_by_group ?: ""
    _itol_outputs                = generate_itol_files(
                                    channel.value(cogs_csv)
                                        .combine(channel.value(color_group_opt))
                                        .combine(synteny_context_data)
                                        .combine(presence_binary_data)
                                        .combine(best_goi_tsv_ch)
                                        .combine(syntenyTaxOrg_script_ch)
                                )


    formatted_fasta_ch          = reformat_fasta_headers(
                                    fasta_ch
                                        .combine(all_genomes_dir)
                                        .combine(new_formatHeaders_script_ch)
                                )
    


    _taxonomy_itol_files_results = taxonomy_itol_files(
                                    taxonomy_info_results.taxonomy_tsv
                                        .combine(df_complete_ch)
                                        .combine(formatted_fasta_ch)
                                        .combine(taxonomy_itol_files_script_ch)
                                )

    
    _maketree_results            = make_tree(formatted_fasta_ch, base_name)
}
