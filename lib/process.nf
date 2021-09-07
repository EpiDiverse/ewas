#!/usr/bin/env nextflow

// parse the samples file, separate the columns, and transpose
process "parsing" {

    label "low"
    label "finish"
    tag "${samples.getName()}"

    maxForks "${params.fork}".toInteger()
    
    publishDir "${params.output}/", pattern: "input/env.txt" , mode: 'copy', enabled: params.input ? true : false
    publishDir "${params.output}/", pattern: "input/cov.txt" , mode: 'copy', enabled: params.input ? true : false
    publishDir "${params.output}/", pattern: "input/gxe.txt" , mode: 'copy', enabled: params.input ? true : false
    
    input:
    path samples

    output:
    path "input/env.txt"
    path "input/cov.txt"
    path "input/gxe.txt"

    when:
    params.input

    script:
    """
    mkdir input

    # remove non-standard chars and sort
    cat ${samples} | tr -d " " | sed 's/\\r\$//' | sort -k1 > input/sorted.tsv

    # generate env file
    echo -e "ID\\tenv" | cat - input/sorted.tsv | cut -f-2 | datamash transpose > input/env.txt

    # generate cov.txt
    paste <(cut -f1 input/sorted.tsv) <(cut -f3- input/sorted.tsv) |
    awk -v OFS="\\t" 'NR==1{printf "ID"; for(i=1; i<=NF-1; i++){printf "%scov%s", OFS,i}; print null}1' |
    datamash transpose > input/cov.txt

    # generate gxe.txt
    cat input/cov.txt <(tail -1 input/env.txt) > input/gxe.txt
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
    
    publishDir "${params.output}/regions", pattern: "${model}/${context}.region.filtered_${model == "Emodel" ? "${params.Emodel_pv}" : model == "Gmodel" ? "${params.Gmodel_pv}" : "${params.GxE_pv}"}_pval.txt.gz" , mode: 'copy', enabled: params.input && (params.DMRs || params.merge) ? true : false
    publishDir "${params.output}/positions", pattern: "${model}/${context}.DMRs.filtered_${model == "Emodel" ? "${params.Emodel_pv}" : model == "Gmodel" ? "${params.Gmodel_pv}" : "${params.GxE_pv}"}_pval.txt.gz" , mode: 'copy', enabled: params.input ? true : false
    publishDir "${params.output}/positions", pattern: "${model}/${context}.DMPs.filtered_${model == "Emodel" ? "${params.Emodel_pv}" : model == "Gmodel" ? "${params.Gmodel_pv}" : "${params.GxE_pv}"}_pval.txt.gz" , mode: 'copy', enabled: params.input ? true : false
    publishDir "${params.output}/positions", pattern: "${model}/${context}.bedGraph.filtered_${model == "Emodel" ? "${params.Emodel_pv}" : model == "Gmodel" ? "${params.Gmodel_pv}" : "${params.GxE_pv}"}_pval.txt.gz" , mode: 'copy', enabled: params.input ? true : false
 
    
    input:
    tuple val(model), val(key), val(context), val(type), path(results), path(logs)
    //tuple model, key, contexts, types, path(results), path(logs)
    
    output:
    tuple val(model), val(key), val(type), path("${model}/*.txt.gz")
    //tuple model, key, val("${types.unique().join("")}"), path("${model}/*.txt")
    //tuple model, key, val("${types.unique().join("")}"), path("input/*.txt"), path("${model}/${key}.filtered_${params.Emodel_pv}_pval.txt")

    when:
    params.input

    script:
    """
    mkdir tmp input ${model}

    total=\$(cat ${logs} | grep "100.00%" | cut -d " " -f3 | tr -d "," | awk 'BEGIN{c=0} {c+=\$0} END{print c}')
    echo -e "${model == "Emodel" ? "ID" : "ID\\tsnp"}\\tbeta\\tstats\\tpvalue\\tFDR" |
    tee input/header.txt ${model}/${key}.txt.gz ${model}/${key}.filtered_${model == "Emodel" ? "${params.Emodel_pv}" : model == "Gmodel" ? "${params.Gmodel_pv}" : "${params.GxE_pv}"}_pval.txt.gz

    # calculate FDR
    if [[ \$(head ${results}.gz | wc -l) == 0 ]]; then
    echo "No findings within current parameter scope" > ${model}/${key}.txt.gz
    else
    sort -T tmp --parallel=${task.cpus} -grk5 ${results}.gz | cut -f${model == "Emodel" ? "2-" : "1-"} |
    awk -F "\\t" -v t="\$total" 'BEGIN{OFS="\\t";p=1;r=t} {fdr=(t/r)*${model == "Emodel" ? "\$4" : "\$5"};
    if(fdr>p){fdr=p}; if(p<=${model == "Emodel" ? "${params.Emodel_pv}" : model == "Gmodel" ? "${params.Gmodel_pv}" : "${params.GxE_pv}"}){print \$0,fdr >> "${model}/${key}.gz.unsorted"};
    print \$0,fdr; p=fdr;r--}' >> ${model}/${key}.txt.gz || exit \$?
    fi

    # sort filtered output
    if [ -f ${model}/${key}.gz.unsorted ]; then
    sort -gk5 ${model}/${key}.gz.unsorted >> ${model}/${key}.filtered_${model == "Emodel" ? "${params.Emodel_pv}" : model == "Gmodel" ? "${params.Gmodel_pv}" : "${params.GxE_pv}"}_pval.txt.gz
    fi
    """ 
}


// GEM_Emodel.out[0]
// process to generate Q-Q plots from Emodel
process "qqPlot" {

    label "low"
    label "ignore"
    tag "${key}"
    
    publishDir "${params.output}/positions", pattern: "${model}/CpG.bedGraph*.png" , mode: 'copy', \
    enabled: params.input && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel) ? true : false
    
    publishDir "${params.output}/positions", pattern: "${model}/CpG.DMRs*.png" , mode: 'copy', \
    enabled: params.input && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel) ? true : false
    
    publishDir "${params.output}/regions", pattern: "${model}/CpG.region*.png" , mode: 'copy', \
    enabled: params.input && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel) ? true : false

    publishDir "${params.output}/positions", pattern: "${model}/CHG.bedGraph*.png" , mode: 'copy', \
    enabled: params.input && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel) ? true : false
    
    publishDir "${params.output}/positions", pattern: "${model}/CHG.DMRs*.png" , mode: 'copy', \
    enabled: params.input && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel) ? true : false
    
    publishDir "${params.output}/regions", pattern: "${model}/CHG.region*.png" , mode: 'copy', \
    enabled: params.input && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel) ? true : false
    
    publishDir "${params.output}/positions", pattern: "${model}/CHH.bedGraph*.png" , mode: 'copy', \
    enabled: params.input && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel) ? true : false
    
    publishDir "${params.output}/positions", pattern: "${model}/CHH.DMRs*.png" , mode: 'copy', \
    enabled: params.input && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel) ? true : false
    
    publishDir "${params.output}/regions", pattern: "${model}/CHH.region*.png" , mode: 'copy', \
    enabled: params.input && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel) ? true : false    
    
    
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
    Rscript ${baseDir}/bin/QQplot.R ${key}.txt.gz ${model}/${key}
    """ 
}
