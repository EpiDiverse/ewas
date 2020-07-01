#!/usr/bin/env nextflow

// filter individual bedGraphs for coverage, filter pairwise comparisons for significance
process "filtering" {

    label "low"
    label "finish"
    tag "${context}.${type}"

    maxForks "${params.fork}".toInteger()
   
    input:
    tuple context, type, sample, path(bed)
    // eg. [CpG, DMPs, sample, /path/to/CpG.bedGraph]

    output:
    tuple context, type, sample, path("${sample}.filtered.bed")
    // eg. [CpG, DMPs, sample, sample.txt]

    when:
    params.input

    script:
    if ( type == "bedGraph" )
        """
        tail -n+2 ${bed} |
        awk 'BEGIN{OFS=\"\\t\"} {if((\$5+\$6)>=${params.coverage}) {printf \"%s\\t%s\\t%s\\t%1.2f\\n\", \$1,\$2,\$3,(\$4/100)}}' |
        sort -k1,1 -k2,2n > ${sample}.filtered.bed
        """  
    else
        """
        file=${workflow.profile.contains("test") ? "${bed}" : "${bed}/${bed}.bed"}
        awk 'BEGIN{OFS=\"\\t\"} \$4<=${params.filter_FDR}{printf \"%s\\t%s\\t%s\\t%1.3f\\n\", \$1,\$2,\$3,\$4}' \$file |
        sort -k1,1 -k2,2n > ${sample}.filtered.bed
        """  
} 


// bedGraph_combined = filtering.out.filter{it[1] == "bedGraph"}.groupTuple()
// DMPs_combined = filtering.out.filter{it[1] == "DMPs"}.groupTuple()
// DMRs_combined = filtering.out.filter{it[1] == "DMRs"}.groupTuple()


//bedtools_unionbedg_input.filter{ it[3].size() > 1 }
//bedtools_unionbedg_input.filter{ it[3].size() == 1 }
// process "process_bedtools_unionbedg_methcalls" { combine filtered bedGraphs
process "bedtools_unionbedg" {

    label "low"
    label "finish"
    tag "${context}.${types.unique().join("")}"

    maxForks "${params.fork}".toInteger()
   
    input:
    tuple context, types, samples, path(beds)
    // eg, [CpG, [DMRs, DMRs, DMRs, ...], [sample1, sample2, sample3, ...], [path1, path2, path3, ...]]

    output:
    tuple context, val("${types.unique().join("")}"), samples, path("${context}.${types.unique().join("")}.bed")
    // eg. [CpG, [DMRs, DMRs, DMRs, ...], [sample1, sample2, sample3, ...], /path/to/CpG.DMRs.bed]

    when:
    params.input

    script:
    """
    bedtools unionbedg -filler NA -i ${beds} -header -names ${samples.join(" ")} > ${context}.${types.unique().join("")}.bed
    """
} 


//bedtools_unionbedg_input.mix(bedtools_unionbedg_output)
// process "process_bedtools_unionbedg_methcalls" { combine filtered bedGraphs
process "bedtools_filtering" {

    label "low"
    label "finish"
    tag "${context}.${type}"

    maxForks "${params.fork}".toInteger()
   
    input:
    tuple context, type, samples, path(bed)
    // eg, [CpG, DMRs, [sample1, sample2, sample3, ...], /path/to/DMRs.bed]

    output:
    tuple context, type, samples, path("bed/${context}.${type}.bed")
    // eg. [CpG, DMRs, [sample1, sample2, sample3, ...], /path/to/DMRs.bed]

    when:
    params.input

    script:
    """
    mkdir bed
    ${samples.getClass() == nextflow.util.ArrayBag && samples.size() > 1 ? "head -1 ${bed}" : "echo -e chrom\\tstart\\tend\\t${samples.join('')}" } \\
    > bed/${context}.${type}.txt

    tail -n+2 ${bed} | awk 'NR!=1{NA=0;c=0;s=0;ss=0;
    for(i=4;i<=NF;i++){if(\$i!="NA"){c++;s+=int(\$i*100+0.5);ss+=int(\$i*100+0.5)^2}else{NA++}};
    sd=sqrt((ss-s^2/c)/c)/100; if(sd>${params.filter_SD} && (NA/(NF-3))<=${params.filter_NA}){print}}' >> bed/${context}.${type}.txt
    ${baseDir}/bin/beta_impute.py bed/${context}.${type}.txt --filter_NA ${filter_NA} --filter_SD ${filter_SD}  > bed/${context}.${type}.bed
    """
} 



//bedtools_unionbedg_input.mix(bedtools_unionbedg_output)
// process "process_bedtools_unionbedg_methcalls" { combine filtered bedGraphs
process "bedtools_sorting" {

    label "low"
    label "finish"
    tag "${context}.${type}"

    maxForks "${params.fork}".toInteger()
   
    input:
    tuple context, type, samples, path(bed)
    // eg, [CpG, DMRs, [sample1, sample2, sample3, ...], /path/to/DMRs.bed]

    output:
    tuple context, type, path("bed/${context}.${type}.bed")
    // eg. [CpG, DMRs, /path/to/DMRs.bed]

    when:
    params.input

    script:
    """
    mkdir tmp bed
    head -1 ${bed} > bed/${context}.${type}.bed
    tail -n+2 ${bed} | sort -T tmp --parallel=${task.cpus} -k1,1 -k2,2n >> bed/${context}.${type}.bed
    """
} 



// bedGraph_DMPs = bedtools_unionbedg.out.filter{it[1] == "bedGraph"}.combine(bedtools_unionbedg.out.filter{it[1] == "DMPs"}, by: 0)
// bedGraph_DMRs = bedtools_unionbedg.out.filter{it[1] == "bedGraph"}.combine(bedtools_unionbedg.out.filter{it[1] == "DMRs"}, by: 0)
// intersect_channel = bedGraph_DMPs.mix(bedGraph_DMRs)
// eg. [CpG, bedGraph, unionbedg.txt, DMPs, unionbedg.txt]


// intersecting bedtools_unionbedg of sample bedGraphs, with DMPs/DMRs
process "bedtools_intersect" {

    label "low"
    label "finish"
    tag "${context}.${type}"

    maxForks "${params.fork}".toInteger()
   
    input:
    tuple context, bedGraph, path("methylation.txt"), type, path("differential.txt")
    // eg. [CpG, bedGraph, CpG.bedGraph.bed, DMPs, CpG.DMPs.bed]

    output:
    tuple context, type, path("${context}.${type}.bed")

    when:
    params.input

    script:
    """
    bedtools intersect -a methylation.txt -b differential.txt -sorted -header > ${context}.${type}.bed
    """  
} 

// bedGraph_DMRs

// takes the same input as previous process bedtools_intersect (bedGraph_DMRs)
// filter regions according to bootstrap support values
process "filter_regions" {

    label "low"
    label "finish"
    tag "${context}.${type}"

    maxForks "${params.fork}".toInteger()
   
    input:
    tuple context, bedGraph, path(methylation), type, path(differential)
    // eg. [CpG, bedGraph, CpG.bedGraph.bed, DMRs, CpG.DMRs.bed]

    output:
    tuple context, bedGraph, path(methylation), val("region"), path("filtered.txt")
     
    when:    
    params.input

    script:
    """
    awk -F "\\t" 'BEGIN{OFS="\\t"} {count=NF-3; for(i=4; i<=NF; i++) {if(\$i=="NA") {count--}};
    if((count/(NF-3))>=${params.proportion}) {print \$0}}' ${differential} > filtered.txt
    """  
} 

// filter_regions.out

// takes the same input as previous process bedtools_intersect (bedGraph_DMRs)
// merging sub-regions identified by bedtools_unionbedg of DMRs
process "bedtools_merge" {

    label "low"
    label "finish"
    tag "${context}.${type}"

    maxForks "${params.fork}".toInteger()
   
    input:
    tuple context, bedGraph, path(methylation), type, path(differential)
    // eg. [CpG, bedGraph, CpG.bedGraph.bed, DMRs, CpG.DMRs.bed]

    output:
    tuple context, bedGraph, path(methylation), val("merged"), path("${context}.merged.txt")
     
    when:    
    params.input && params.merge

    script:
    """
    cat <(head -1 ${differential} | cut -f-3) <(bedtools merge -i ${differential} | sort -k1,1 -k2,2n) > ${context}.merged.txt
    """

} 


// average_channel = filter_regions.out.mix(bedtools_merge.out)
// eg. [CpG, bedGraph, methylation.txt, DMPs, differential.txt] 

// takes the same input as previous process bedtools_intersect (bedGraph_DMRs)
// intersecting bedtools_unionbedg of sample bedGraphs, with 
process "average_over_regions" {

    label "low"
    label "finish"
    tag "${context}.${type}"

    maxForks "${params.fork}".toInteger()
   
    input:
    tuple context, bedGraph, path("methylation.txt"), type, path("differential.txt")
    // eg. [CpG, bedGraph, CpG.bedGraph.bed, DMRs, CpG.DMRs.bed]

    output:
    tuple context, type, path("${context}.${type}.bed")

    when:
    params.input

    script:
    """
    tail -q -n+2 differential.txt methylation.txt | cut -f1 | uniq | sort | uniq > index.txt
    ${baseDir}/bin/average_over_bed.py <(tail -n+2 differential.txt) methylation.txt index.txt > ${context}.${type}.bed
    """

} 


// meth_channel = bedtools_unionbedg.out.filter{it[1] == "bedGraph"}.mix(bedtools_intersect.out, average_over_regions.out)
// eg. [CpG, bedGraph, methylation.txt]
// eg. [CpG, DMRs, CpG.bed]
// eg. [CpG, regions, CpG.avg]


// RUN GEM Emodel
process "GEM_Emodel" {
    
    label "low"
    label "finish"
    tag "${context}.${type} - ${meth.baseName}"

    input:
    tuple context, type, path(meth)
    path envs
    path covs
    
    output:
    //tuple context, type, path("output/*.txt"), path("output/*.log")
    path "output/${context}.${type}.gz"
    path "output/${context}.${type}.log"
   
    when:
    params.input && (!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel
    
    script: 
    """
    mkdir output
    awk -F "\\t" '{printf \"%s:%s-%s\",\$1,\$2,\$3; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' ${meth} > \$(basename ${meth} .bed).txt
    Rscript ${baseDir}/bin/GEM_Emodel.R ${baseDir}/bin ${envs} ${covs} \$(basename ${meth} .bed).txt ${params.Gmodel_pv} output/temp > output/${context}.${type}.log || exit \$?
    tail -n+2 output/temp.txt | gzip > output/${context}.${type}.gz && rm output/temp.txt
    """

}



// GEM_Emodel.out[0]
// process to generate manhattan plots from Emodel
process "manhattan" {

    label "low"
    label "ignore"
    tag "${key}"
     
    input:
    tuple model, key, type, path(txt)
    // eg. [Emodel, CpG.bedGraph, bedGraph, [/paths/... ,/paths/...]]
    
    output:
    tuple type, path("*.png") optional true
    tuple type, path("*.zip") optional true

    when:
    params.input && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel)

    script:
    """
    awk -F "\\t" 'BEGIN{OFS="\\t"; print "SNP","CHR","BP","P"} NR!=1{split(\$1,cpg,":"); split(cpg[2],cpos,"-"); pos=(cpos[1]+cpos[2])/2;
    print \$1,cpg[1],pos,\$5}' ${key}.filtered_${params.output_FDR}_FDR.txt > manhattan.txt
    Rscript ${baseDir}/bin/manhattan.R manhattan.txt ${key}.filtered_${params.output_FDR}_FDR 0.00000001 0.000001
    """ 
}
