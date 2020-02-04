#!/usr/bin/env nextflow

//initial vars

if (params.DMP && params.DMR) {
    println 'ERROR: Please specify at most one input type!'
}
else if (params.DMP) {   
    CpG_metilene_path        = "${params.DMP}/CpG/metilene/*" 
    CHG_metilene_path        = "${params.DMP}/CHG/metilene/*"
    CHH_metilene_path        = "${params.DMP}/CHH/metilene/*"   
}
else if (params.DMR) {   
    CpG_metilene_path        = "${params.DMR}/CpG/metilene/*" 
    CHG_metilene_path        = "${params.DMR}/CHG/metilene/*"
    CHH_metilene_path        = "${params.DMR}/CHH/metilene/*"   
}


CpG_path = "${params.input}/{${commaLine[0..-2]}}/bedGraph/*CpG.bedGraph"
CHG_path = "${params.input}/{${commaLine[0..-2]}}/bedGraph/*CHG.bedGraph"
CHH_path = "${params.input}/{${commaLine[0..-2]}}/bedGraph/*CHH.bedGraph"


//CpG_path                 = ["${params.input}/CpG/metilene/*.bed", "{params.input}/CpG/input/*.txt"]
//CpG_path                 = tuple("${params.input}/CpG/metilene/*[\\d].[\\d+].bed", "${params.input}/CpG/input/*.txt")
//CpG_path                 = tuple("${params.input}/CpG/metilene/*[0-9].[0-9][0-9].bed", "${params.input}/CpG/input/*.txt")



// STAGE DMP,DMR and METHCALLS INPUT DIRECTORIES 


//METHCALLS CHANNELS FOR ALL CONTEXTS
//eg [group1_vs_group2.0.05.bed]
bedGraph_channel1 = params.noCpG  ? Channel.empty() :
Channel
    .fromFilePairs( CpG_path, size: 1)
    .ifEmpty{ exit 1, "ERROR: 2: ${params.input}\n"}
    .set{CpG_input1}

//bedGraph_channel1.into{CpG_input1; CpG_input2; CpG_input3}
methcalls1 = (params.noCpG  ? Channel.empty() :CpG_input1.combine(samples_bedGraph, by: 0).map{tuple("CpG", *it)})

//CHG input channel
bedGraph_channel2 = params.noCHG  ? Channel.empty() :
Channel
    .fromFilePairs( CHG_path, size: 1)
    .ifEmpty{ exit 1, "ERROR: 2: ${params.input}\n"}
    .set{CHG_input1}
 
methcalls2 =(params.noCHG  ? Channel.empty() : CHG_input1.combine(samples_bedGraph2, by: 0).map{tuple("CHG", *it)})

//CHH input channel
bedGraph_channel3 = params.noCHH  ? Channel.empty() :
Channel
    .fromFilePairs( CHH_path, size: 1)
    .ifEmpty{ exit 1, "ERROR: 2: ${params.input}\n"}
    .set{CHH_input1}
methcalls3 = (params.noCHH  ? Channel.empty() :CHH_input1.combine(samples_bedGraph3, by: 0).map{tuple("CHH", *it)})

methcalls= methcalls1.mix(methcalls2).mix(methcalls3)
//methcalls.view()


//DMP CHANNELS FOR ALL CONTEXTS
//CpG DMP channel
metilene_channel1= params.noCpG || !params.DMP  ? Channel.empty() :
Channel
    .fromFilePairs( CpG_metilene_path,  type: 'dir', size:1)
    .ifEmpty{ exit 1, "ERROR: No input found in: ${params.DMP}\n" }

DMPs1 = (params.noCpG  ? Channel.empty() : metilene_channel1.map{tuple("CpG", *it)}.groupTuple())

metilene_channel3= params.noCHG  || !params.DMP ? Channel.empty() :
Channel
    .fromFilePairs( CHG_metilene_path, type: 'dir', size: 1)
    .ifEmpty{ exit 1, "ERROR: 1: ${params.DMP}\n" }
DMPs2 = (params.noCHG  ? Channel.empty() : metilene_channel3.map{tuple("CHG", *it)}.groupTuple())


metilene_channel5= params.noCHH  || !params.DMP ? Channel.empty() :
Channel
    .fromFilePairs( CHH_metilene_path, type: 'dir', size: 1)
    .ifEmpty{ exit 1, "ERROR: 1: ${params.DMP}\n" }
DMPs3 = (params.noCHH  ? Channel.empty() :metilene_channel5.map{tuple("CHH", *it)}.groupTuple())

DMPs= DMPs1.mix(DMPs2).mix(DMPs3)



//DMR CHANNELS FOR ALL CONTEXTS
//CpG DMR channel
metilene_channel2= params.noCpG  || !params.DMR ? Channel.empty() :
Channel
    .fromFilePairs( CpG_metilene_path, type: 'dir', size: 1)
    .ifEmpty{ exit 1, "ERROR: No input found in: ${params.DMR}\n" }

DMRs1 = (params.noCpG  ? Channel.empty() :metilene_channel2.map{tuple("CpG", *it)}.groupTuple())


metilene_channel4= params.noCHG  || !params.DMR ? Channel.empty() :
Channel
    .fromFilePairs( CHG_metilene_path, type: 'dir', size: 1)
    .ifEmpty{ exit 1, "ERROR: 1: ${params.DMR}\n" }
DMRs2 = (params.noCHG  ? Channel.empty() : metilene_channel4.map{tuple("CHG", *it)}.groupTuple())


metilene_channel6= params.noCHH  || !params.DMR ? Channel.empty() :
Channel
    .fromFilePairs( CHH_metilene_path, type: 'dir', size: 1)
    .ifEmpty{ exit 1, "ERROR: 1: ${params.DMR}\n" }
DMRs3 =(params.noCHH  ? Channel.empty() : metilene_channel6.map{tuple("CHH", *it)}.groupTuple())

DMRs= DMRs1.mix(DMRs2).mix(DMRs3)



// BEGIN PIPELINE

// taking sample.tsv and generating env.txt and cov.txt for all input types (meth. calls, DMPs and DMRs)

process "process_cov_and_env_files" {

    label "low"
    publishDir "${output_path}/meth_calls", mode: 'copy'
    if(params.DMP){publishDir "${output_path}/DMPs", mode: 'copy'}
    if(params.DMR){publishDir "${output_path}/DMRs", mode: 'copy'}

    input: 
    file samples from samples_file

    output:   
    file "env.txt" into environment_methcalls, environment_DMPs, environment_DMRs, env,gxe
    file "cov.txt" into covariate_methcalls, covariate_DMPs, covariate_DMRs, cov_gxe
                  
    script:
    """
    #!/bin/bash
    sed '/[0-9]\\,/s/\\,/./g' ${samples} | awk -F "\\t" '{printf \"%s\t%.6f",\$1,\$2; for(i=3; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' | awk -F "\\t" '{printf \"%s\\t%s\",\$1,\$2; for(i=3; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' > samples2.txt
    cat samples2.txt |grep '[0-9]'| sed 's/,//g' > samples3.txt
    echo -e "ID\tenv" | cat - samples3.txt > samples4.txt
    cat <(cut -f1 samples4.txt | paste -s) <(cut -f2 samples4.txt | paste -s) > env.txt
    cut -f1 samples3.txt > header.txt
    cut -d \$'\\t' -f3- samples3.txt  > pre_cov.txt
    paste header.txt pre_cov.txt > 2pre_cov.txt
     awk 'NR==1{printf "ID"; for(i=1; i<=NF-1; i++) h=h OFS "cov" i; print h}1' OFS='\t' 2pre_cov.txt > 3pre_cov.txt
    cat 3pre_cov.txt | datamash transpose | tr -d "\\r" > cov.txt    
   
   """
}

/*

process "process_gxe_file" {

    label "low"
    publishDir "${output_path}/GEM_GXEmodel/input", mode: 'copy'
   
    input:
    file cov from cov_gxe
    file env from env_gxe

    output:
    file("gxe.txt") into gxe_methcalls, gxe_DMPs, gxe_DMRs
    
    when:
    params.GEM_GXEmodel
    
    script:
    """
    #!/bin/bash
    head -n 1 ${env} > 2env.txt
    paste ${cov} 2env.txt > gxe.txt
    """  
} 

*/

process "process_input_files_meth_calls" {

    label "low"
    //publishDir "${output_path}/meth_calls/${context}/GEM_Emodel/input", mode: 'copy'
   
    input:
    set context, sample, file(bedGraph), env, cov from methcalls

    output:
    set context, sample, file("*.txt"), env, cov into cut1, cut2, cut3
    // eg. [CpG, samplename, bedgraph.txt]
    

    script:
    """
    #!/bin/bash
    tail -n+2 ${bedGraph} | awk 'BEGIN{OFS=\"\\t\"} {if((\$5+\$6)>=${params.coverage}) {printf \"%s\\t%s\\t%s\\t%1.2f\\n\", \$1,\$2,\$3,(\$4/100)}}' > ${bedGraph}.txt
    """  
} 


//DMP INPUT SECTION
// taking lists of files with with DMPs and merging into GEM E_model methylation input format


process "process_bedtools_unionbedg_methcalls1" {

    label "low"
    publishDir "${output_path}/DMPs/${context}/GEM_Emodel/input", mode: 'copy'
   
    input:   
    set context, samples, file(bedGraph), env, cov  from cut1.groupTuple()

    output:
    set context, samples, file("methylation.txt"), env, cov into bedGraph_DMPs
    
    when:
    params.DMP

    script:
    """
    #!/bin/bash
    bedtools unionbedg -filler NA -i ${bedGraph} -header -names ${samples.join(" ")} | awk -F "\\t" '{printf \"%s\t%s\t%s",\$1,\$2,\$3; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' | sed 's/NA/ /g' > methylation.txt 
 
    """  
} 


process "process_bedtools_unionbedg_methcalls2" {

    label "low"
    //publishDir "${output_path}/meth_calls/${context}/GEM_Emodel/input", mode: 'copy'
   
    input:
   
    set context, samples, file(bedGraph), env, cov  from cut2.groupTuple()

    output:
    set context, samples, file("methylation2.txt"), env, cov into bedGraph_methcalls
    
    when:
    !params.DMP && !params.DMR
        
    script:
    """
    #!/bin/bash
    bedtools unionbedg -filler NA -i ${bedGraph} -header -names ${samples.join(" ")} | awk -F "\\t" '{printf \"%s_%s\",\$1,\$2; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' | sed 's/NA/ /g' > methylation2.txt 
    
    """  
} 

process "process_bedtools_unionbedg_methcalls3" {

    label "low"
    //publishDir "${output_path}/DMRs/${context}/GEM_Emodel/input", mode: 'copy'
   
    input:   
    set context, samples, file(bedGraph), env, cov  from cut3.groupTuple()

    output:
    set context, samples, file("methylation.txt"), env, cov into bedGraph_DMRs
    
    when:
    params.DMR

    script:
    """
    #!/bin/bash
    bedtools unionbedg -filler NA -i ${bedGraph} -header -names ${samples.join(" ")} | awk -F "\\t" '{printf \"%s\t%s\t%s",\$1,\$2,\$3; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' | sed 's/NA/ /g' > methylation.txt 
    
    """  
} 



//use this inside GEM for methcalls
//bedtools unionbedg -filler NA -i ${bed} -header -names ${samples.join(" ")} | awk -F "\\t" '{printf \"%s_%s\",\$1,\$2; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' | sed 's/NA/ /g' > methylation2.txt 



// RUN GEM with methylation calls
process "process_GEM_Emodel_run_meth_calls" {
    
    label "low"    
    publishDir "${output_path}/meth_calls/${context}/GEM_Emodel", mode: 'copy'
    
    //beforeScript "${workflow.profile == "standard" || workflow.profile == "diverse" ? 'source activate /scr/epi/pipelines/ewas/libs/gem' : ''}"
    //afterScript "${workflow.profile == "standard" || workflow.profile == "diverse" ? 'conda deactivate /scr/epi/pipelines/ewas/libs/gem' : ''}"

    input:
    file envs from environment_methcalls
    file covs from covariate_methcalls
    set context, samples, file (meth), env, cov from bedGraph_methcalls
    
    output:
    file("emodel.txt")
    file("filtered_emodel.txt")
    file("emodel.png") 

    when:    
    !params.DMP && !params.DMR
    
    script: 
    """
    Rscript ${baseDir}/bin/GEM_Emodel_methcalls.R ${envs} ${covs} ${meth} ${params.Emodel_pv} > emodel.txt
    sort -n emodel.txt | awk 'BEGIN{OFS=\"\\t\"} {if((\$4)<=${params.sig} && (\$5)<=${params.FDR}) {printf \"%s\\t%s\\t%s\\t%s\\t%s\\n\", \$1,\$2,\$3,\$4,\$5}}' > filtered_emodel.txt    
    """
}

/*

//
//
//
//
//
//DMP PROCESS STARTS
process "process_bedtools_unionbedg_DMPs" {

    label "low"
    //publishDir "${output_path}/DMPs/${context}/GEM_Emodel/input", mode: 'copy'
   
    input:
    set context, pairwise, file(bed) from DMPs.map{tuple(it[0],it[1],tuple(it[2].flatten()))}
   
    output:
    
    set context, pairwise, file("dmp.txt")  into pre_DMP 
    
    when:    
    params.DMP
   
    script:
    """
    #!/bin/bash
    
    find . -mindepth 1 -maxdepth 1 -type l | while read dir; do file=\$(ls \$dir/*.bed | awk 'BEGIN{OFS="\\t"} {print length,\$0}' | sort -nrk1 | head -1 | cut -f2); cut -f1,2,3,4 \$file > \$(basename \$file); done
    bedtools unionbedg -i *.bed > 2pre_dmp.txt
    cut -f1,2,3,4 2pre_dmp.txt > dmp.txt
    
    """  
} 



process "process_bedtools_intersect_DMP" {

    label "low"
    //publishDir "${output_path}/DMPs/${context}/GEM_Emodel/input", mode: 'copy'
   
    input:
      
    set context, pairwise, file (dmp) from pre_DMP
    set context, samples, file (meth), env, cov from bedGraph_DMPs
   
    output:
    set context, pairwise, file("methylation_DMP.txt") into methylation_DMP 
     
    when:    
    params.DMP
   
    script:
    """
    #!/bin/bash
    head -1 ${meth} > header.txt
    sed '1d' ${meth} > 2meth.txt
    awk '{if(\$2<\$3) {printf \"\\t%s\",\$i}; print null}' ${dmp} > 2dmp.txt
    bedtools intersect -a 2meth.txt -b 2dmp.txt > pre_dmp.txt
    cat header.txt pre_dmp.txt > 2pre_dmp.txt 
    cat 2pre_dmp.txt | awk -F "\\t" '{printf \"%s_%s\",\$1,\$2; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' > methylation_DMP.txt
 
    """  
} 
//cat header.txt 2meth.txt > methylation_DMP.txt
//    bedtools intersect -a 2meth.txt -b ${dmp} > 3meth.txt
//    cat header.txt 3meth.txt > methylation_DMP.txt

//file("methylation.txt") from bedGraph_DMPs
//   bedtools intersect -a methylation.txt -b dmp.txt > pre2meth.txt
//    cat header.txt pre2meth.txt > methylation_DMP.txt



process "process_GEM_Emodel_run_DMPs" {
    
    label "low"    
    publishDir "${output_path}/DMPs/${context}/GEM_Emodel", mode: 'copy'
    

    input:
    file envs from environment_DMPs
    file covs from covariate_DMPs
    set context, pairwise, file (meth) from methylation_DMP

    output:
    file("emodel.txt")
    file("filtered_emodel.txt")
    file("emodel.png") 

    when: 
    params.DMP
    
    script: 
    """
    Rscript ${baseDir}/bin/GEM_Emodel_methcalls.R ${envs} ${covs} ${meth} ${params.Emodel_pv} > emodel.txt
    sort -n emodel.txt | awk 'BEGIN{OFS=\"\\t\"} {if((\$4)<=${params.sig} && (\$5)<=${params.FDR}) {printf \"%s\\t%s\\t%s\\t%s\\t%s\\n\", \$1,\$2,\$3,\$4,\$5}}' > filtered_emodel.txt        
    """
}

*/

//
//
//
//
//
//DMR PROCESS STARTS
process "process_bedtools_unionbedg_DMRs" {

    label "low"
    //publishDir "${output_path}/DMRs/${context}/GEM_Emodel/input", mode: 'copy'
   
    input:
    set context, pairwise, file(bed) from DMRs.map{tuple(it[0],it[1],tuple(it[2].flatten()))}
   
    output:
    
    set context, pairwise, file("dmr.txt")  into pre_DMR 
    
    when:    
    params.DMR
   
    script:
    """
    #!/bin/bash
    
    find . -mindepth 1 -maxdepth 1 -type l | while read dir; do file=\$(ls \$dir/*.bed | awk 'BEGIN{OFS="\\t"} {print length,\$0}' | sort -nrk1 | head -1 | cut -f2); cut -f1,2,3,4 \$file > \$(basename \$file); done
    bedtools unionbedg -i *.bed > 2pre_dmr.txt
    cut -f1,2,3,4 2pre_dmr.txt > dmr.txt
    
    """  
} 



process "process_bedtools_intersect_DMR" {

    label "low"
    //publishDir "${output_path}/DMRs/${context}/GEM_Emodel/input", mode: 'copy'
   
    input:
      
    set context, pairwise, file (dmr) from pre_DMR
    set context, samples, file (meth), env, cov from bedGraph_DMRs
   
    output:
    set context, pairwise, file("methylation_DMR.txt") into methylation_DMR 
     
    when:    
    params.DMR
   
    script:
    """
    #!/bin/bash
    head -1 ${meth} > header.txt
    sed '1d' ${meth} > 2meth.txt
    awk '{if(\$2<\$3) {print \$0}}' ${dmr} > 2dmr.txt
    bedtools intersect -a 2meth.txt -b 2dmr.txt > pre_dmr.txt
    cat header.txt pre_dmr.txt > 2pre_dmr.txt 
    cat 2pre_dmr.txt | awk -F "\\t" '{printf \"%s_%s\",\$1,\$2; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' > methylation_DMR.txt
 
    """  
} 

process "process_GEM_Emodel_run_DMRs" {
    
    label "low"    
    publishDir "${output_path}/DMRs/${context}/GEM_Emodel", mode: 'copy'
    

    input:
    file envs from environment_DMRs
    file covs from covariate_DMRs
    set context, pairwise, file (meth) from methylation_DMR

    output:
    file("emodel.txt") 
    file("emodel.png") 
    file("filtered_emodel.txt")

    when:   
    params.DMR
    
    script: 
    """
    Rscript ${baseDir}/bin/GEM_Emodel_methcalls.R ${envs} ${covs} ${meth} ${params.Emodel_pv} > emodel.txt
    sort -n emodel.txt | awk 'BEGIN{OFS=\"\\t\"} {if((\$4)<=${params.sig} && (\$5)<=${params.FDR}) {printf \"%s\\t%s\\t%s\\t%s\\t%s\\n\", \$1,\$2,\$3,\$4,\$5}}' > filtered_emodel.txt    
   
    """
}
