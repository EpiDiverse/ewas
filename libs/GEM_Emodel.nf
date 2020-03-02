#!/usr/bin/env nextflow

// filter individual bedGraphs for coverage, filter pairwise comparisons for significance
process "filtering" {

    label "low"
    label "finish"
    tag "$type - $context"

    maxForks "${params.fork}".toInteger()
   
    input:
    tuple context, type, sample, path(bed)
    // eg. [CpG, DMPs, sample, /path/to/CpG.bedGraph]

    output:
    tuple context, type, sample, path("${sample}.bed")
    // eg. [CpG, DMPs, sample, sample.txt]

    when:
    params.input

    script:
    if ( type == "bedGraph" )
        """
        tail -n+2 ${bed} |
        awk 'BEGIN{OFS=\"\\t\"} {if((\$5+\$6)>=${params.coverage}) {printf \"%s\\t%s\\t%s\\t%1.2f\\n\", \$1,\$2,\$3,(\$4/100)}}' |
        sort -k1,1 -k2,2n > ${sample}.bed
        """  
    else
        """
        awk 'BEGIN{OFS=\"\\t\"} \$4<=${params.sig}{printf \"%s\\t%s\\t%s\\t%1.3f\\n\", \$1,\$2,\$3,\$4}' ${bed}/${bed}.bed |
        sort -k1,1 -k2,2n > ${sample}.bed
        """  
} 


// bedGraph_combined = filtering.out.filter{it[1] == "bedGraph"}.groupTuple()
// DMPs_combined = filtering.out.filter{it[1] == "DMPs"}.groupTuple()
// DMRs_combined = filtering.out.filter{it[1] == "DMRs"}.groupTuple()


// process "process_bedtools_unionbedg_methcalls" { combine filtered bedGraphs
process "bedtools_unionbedg" {

    label "low"
    label "finish"
    tag "${types.unique().join("")} - $context"

    maxForks "${params.fork}".toInteger()
   
    input:
    tuple context, types, samples, path(beds)
    // eg, [CpG, [DMRs, DMRs, DMRs, ...], [sample1, sample2, sample3, ...], [path1, path2, path3, ...]]

    output:
    tuple context, val("${types.unique().join("")}"), path("${context}.${types.unique().join("")}.bed")
    // eg. [CpG, DMRs, /path/to/DMRs.bed]

    when:
    params.input

    script:
    """
    bedtools unionbedg -filler NA -i ${beds} -header -names ${samples.join(" ")} > unsorted.${context}.${types.unique().join("")}.bed || exit \$?
    head -1 unsorted.${context}.${types.unique().join("")}.bed > ${context}.${types.unique().join("")}.bed
    tail -n+2 unsorted.${context}.${types.unique().join("")}.bed | sort -k1,1 -k2,2n >> ${context}.${types.unique().join("")}.bed
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
    tag "$type - $context"

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
    tag "$type - $context"

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
    if((count/(NF-3))>=${params.bootstrap}) {print \$0}}' ${differential} > filtered.txt
    """  
} 

// filter_regions.out

// takes the same input as previous process bedtools_intersect (bedGraph_DMRs)
// merging sub-regions identified by bedtools_unionbedg of DMRs
process "bedtools_merge" {

    label "low"
    label "finish"
    tag "$type - $context"

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
    bedtools merge -i ${differential} | sort -k1,1 -k2,2n > ${context}.merged.txt
    """  
} 


// average_channel = filter_regions.out.mix(bedtools_merge.out)
// eg. [CpG, bedGraph, methylation.txt, DMPs, differential.txt] 

// takes the same input as previous process bedtools_intersect (bedGraph_DMRs)
// intersecting bedtools_unionbedg of sample bedGraphs, with 
process "average_over_regions" {

    label "low"
    label "finish"
    tag "$type - $context"

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
    ${baseDir}/bin/average_over_bed.py <(tail -n+2 differential.txt) <(cut -f1,3- methylation.txt) > ${context}.${type}.bed
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
    tag "$type - $context"

    input:
    tuple context, type, path(meth)
    path envs
    path covs
    
    output:
    tuple type, path("${context}.${type}.filtered_${params.FDR}_FDR.txt")
    tuple type, path("${context}.${type}.txt")
    tuple type, path("${context}.${type}.jpg")
    tuple type, path("${context}.${type}.log")
   
    when:
    (!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel
    
    script: 
    """
    awk -F "\\t" '{printf \"%s:%s-%s\",\$1,\$2,\$3; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' ${meth} > ${context}.txt
    Rscript ${baseDir}/bin/GEM_Emodel.R ${envs} ${covs} ${context}.txt ${params.Emodel_pv} ${context}.${type} > ${context}.${type}.log
    sort -V ${context}.${type}.txt |
    awk 'BEGIN{OFS="\\t"; print "cpg\\tbeta\\tstats\\tpvalue\\tFDR"} \$5<=${params.FDR}{print \$1,\$2,\$3,\$4,\$5}' > ${context}.${type}.filtered_${params.FDR}_FDR.txt
    """
}