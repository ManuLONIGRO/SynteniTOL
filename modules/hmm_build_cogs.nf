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