#!/usr/bin/env nextflow

// parse the samples file, separate the columns, and transpose
process "parsing" {

    label "low"
    label "finish"
    tag "${samples.getName()}"

    maxForks "${params.fork}".toInteger()
    
    publishDir "${params.output}/input", pattern: "env.txt" , mode: 'copy', enabled: params.input ? true : false
    publishDir "${params.output}/input", pattern: "cov.txt" , mode: 'copy', enabled: params.input ? true : false
    publishDir "${params.output}/input", pattern: "gxe.txt" , mode: 'copy', enabled: params.input ? true : false
    
    input:
    path samples

    output:
    path "env.txt"
    path "cov.txt"
    path "gxe.txt"

    when:
    params.input

    script:
    """
    sed '/[0-9]\\,/s/\\,/./g' ${samples} |
    awk -F "\\t" '{printf \"%s\\t%.6f",\$1,\$2; for(i=3; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' |
    awk -F "\\t" '{printf \"%s\\t%s\",\$1,\$2; for(i=3; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' > samples2.txt
    
    cat samples2.txt | grep '[0-9]'| sed 's/,//g' > samples3.txt
    echo -e "ID\\tenv" | cat - samples3.txt > samples4.txt
    cat <(cut -f1 samples4.txt | paste -s) <(cut -f2 samples4.txt | paste -s) > env.txt

    cut -f1 samples3.txt > header.txt
    cut -d \$'\\t' -f3- samples3.txt  > pre_cov.txt
    
    paste header.txt pre_cov.txt > 2pre_cov.txt
    awk 'NR==1{printf "ID"; for(i=1; i<=NF-1; i++) h=h OFS "cov" i; print h}1' OFS='\\t' 2pre_cov.txt > 3pre_cov.txt
    cat 3pre_cov.txt | datamash transpose | tr -d "\\r" > cov.txt

    cat cov.txt <(tail -1 env.txt) > gxe.txt
    """
}



// meth_channel = bedtools_unionbedg.out.filter{it[1] == "bedGraph"}.mix(bedtools_intersect.out, average_over_regions.out)
// eg. [CpG, bedGraph, methylation.txt]
// eg. [CpG, DMRs, CpG.bed]
// eg. [CpG, regions, CpG.avg]


// bcftools.out
// extract required format
process "split_scaffolds" {

    label "low"
    label "finish"
    tag "${context}.${type}"
     
    input:
    tuple val(context), val(type), path(bed)
    
    output:
    tuple val(context), val(type), path("output/*.bed")

    when:
    params.input

    script:
    """   
    mkdir output
    awk -F "\\t" '{if(NR-1){if(\$1 in arr == 0){arr[\$1]=\$1; print header > "output/"\$1".bed"};
    print \$0 >> "output/"\$1".bed"} else {header=\$0}}' ${bed}
    """ 
}
//split_scaffolds.out.transpose()



// GEM_Gmodel.out.map{ tuple( it[0] + "." + it[1], *it) }.groupTuple().map{ tuple("Gmodel", *it) }
// process to calculate FDR on combined files after splitting
process "calculate_FDR" {

    label "high"
    label "finish"
    tag "${model}:${key}"
    
    publishDir "${params.output}/regions", pattern: "${model}/${context}.region.filtered_${params.output_FDR}_FDR.txt" , mode: 'copy', enabled: params.input && (params.DMRs || params.merge) ? true : false
    publishDir "${params.output}/regions", pattern: "${model}/${context}.region.txt" , mode: 'copy', enabled: params.input && (params.DMRs || params.merge) ? true : false 
    publishDir "${params.output}/positions", pattern: "${model}/${context}.DMRs.filtered_${params.output_FDR}_FDR.txt" , mode: 'copy', enabled: params.input ? true : false
    publishDir "${params.output}/positions", pattern: "${model}/${context}.DMRs.txt" , mode: 'copy', enabled: params.input  ? true : false
    publishDir "${params.output}/positions", pattern: "${model}/${context}.DMPs.filtered_${params.output_FDR}_FDR.txt" , mode: 'copy', enabled: params.input ? true : false
    publishDir "${params.output}/positions", pattern: "${model}/${context}.DMPs.txt" , mode: 'copy', enabled: params.input  ? true : false
    publishDir "${params.output}/positions", pattern: "${model}/${context}.bedGraph.filtered_${params.output_FDR}_FDR.txt" , mode: 'copy', enabled: params.input ? true : false
    publishDir "${params.output}/positions", pattern: "${model}/${context}.bedGraph.txt" , mode: 'copy', enabled: params.input  ? true : false
 
    
    input:
    tuple val(model), val(key), val(context), val(type), path(results), path(logs)
    //tuple model, key, contexts, types, path(results), path(logs)
    
    output:
    tuple val(model), val(key), val(type), path("${model}/*.txt")
    //tuple model, key, val("${types.unique().join("")}"), path("${model}/*.txt")
    //tuple model, key, val("${types.unique().join("")}"), path("input/*.txt"), path("${model}/${key}.filtered_${params.output_FDR}_FDR.txt")

    when:
    params.input

    script:
    """
    mkdir tmp input ${model}
    #tail -q -n+2 ${results} > input/${key}.txt
    total=\$(cat ${logs} | grep "100.00%" | cut -d " " -f3 | tr -d "," | awk 'BEGIN{c=0} {c+=\$0} END{print c}')
    echo -e "${model == "Emodel" ? "cpg" : "cpg\\tsnp"}\\tbeta\\tstats\\tpvalue\\tFDR" |
    tee input/header.txt ${model}/${key}.txt ${model}/${key}.filtered_${params.output_FDR}_FDR.txt

    #if [[ \$(head input/${key}.txt | wc -l) == 0 ]]; then
    if [[ \$(head ${results} | wc -l) == 0 ]]; then
    echo "No findings with ${model == "Emodel" ? "--Emodel_pv ${params.Emodel_pv}" : model == "Gmodel" ? "--Gmodel_pv ${params.Gmodel_pv}" : "--GxE_pv ${params.GxE_pv}"}" > ${model}/${key}.txt
    else
    sort -T tmp --parallel=${task.cpus} -grk5 ${results} | cut -f${model == "Emodel" ? "2-" : "1-"} |
    awk -F "\\t" -v t="\$total" 'BEGIN{OFS="\\t";p=1;r=t} {fdr=(t/r)*${model == "Emodel" ? "\$4" : "\$5"};
    if(fdr>p){fdr=p}; if(fdr<=${params.output_FDR}){print \$0,fdr >> "${model}/${key}.filtered_${params.output_FDR}_FDR.txt"};
    print \$0,fdr; p=fdr;r--}' >> ${model}/${key}.txt || exit \$?
    fi
    """ 
}


// GEM_Emodel.out[0]
// process to generate manhattan plots from Emodel
process "qqPlot" {

    label "low"
    label "ignore"
    tag "${key}"
    
    publishDir "${params.output}/positions", pattern: "${model}/${context}.DMRs*.png" , mode: 'copy', enabled: params.input ? true : false
    publishDir "${params.output}/regions", pattern: "${model}/${context}.region*.png" , mode: 'copy', enabled: params.input && params.DMRs ? true : false
    
    input:
    //tuple val(model), val(key), val(context), val(type), path(result)
    tuple val(model), val(key), val(type), path(result)
    // eg. [Emodel, CpG.bedGraph, bedGraph, [/paths/... ,/paths/...]]
    
    output:
    //tuple val(model), val(key), val(type), path("${model}/*.png")
    tuple val(model), val(key), val(type), path("${model}/*.png") optional true

    when:
    params.input

    script:
    """
    mkdir ${model}
    Rscript ${baseDir}/bin/QQplot.R ${key}.txt ${model}/${key}
    """ 
}
