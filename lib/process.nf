#!/usr/bin/env nextflow

// parse the samples file, separate the columns, and transpose
process "parsing" {

    label "low"
    label "finish"
    tag "${samples.getName()}"

    maxForks "${params.fork}".toInteger()

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
    tuple context, type, path(bed)
    
    output:
    tuple context, type, path("output/*.bed")

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
     
    input:
    tuple model, key, context, type, path(results), path(logs)
    //tuple model, key, contexts, types, path(results), path(logs)
    
    output:
    tuple model, key, type, path("${model}/*.txt")
    //tuple model, key, val("${types.unique().join("")}"), path("${model}/*.txt")
    //tuple model, key, val("${types.unique().join("")}"), path("input/*.txt"), path("${model}/${key}.filtered_${params.output_FDR}_FDR.txt")

    when:
    params.input

    script:
    """
    mkdir tmp input ${model}
    total=\$(cat ${logs} | grep "100.00%" | cut -d " " -f3 | tr -d "," | awk 'BEGIN{c=0} {c+=\$0} END{print c}')
    echo -e "${model == "Emodel" ? "cpg" : "cpg\\tsnp"}\\tbeta\\tstats\\tpvalue\\tFDR" |
    tee input/header.txt ${model}/${key}.txt ${model}/${key}.filtered_${params.output_FDR}_FDR.txt

    if [ -z \$(gzip -cd ${results} | head -c1) ]; then
    echo "No findings with ${model == "Emodel" ? "--Emodel_pv ${params.Emodel_pv}" : model == "Gmodel" ? "--Gmodel_pv ${params.Gmodel_pv}" : "--GxE_pv ${params.GxE_pv}"}" > ${model}/${key}.txt
    else
    gzip -cd ${results} | sort -T tmp --parallel=${task.cpus} -grk5 | cut -f${model == "Emodel" ? "2-" : "1-"} |
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
     
    input:
    tuple model, key, type, path(result)
    // eg. [Emodel, CpG.bedGraph, bedGraph, [/paths/... ,/paths/...]]
    
    output:
    tuple type, path("${model}/*.png") optional true

    when:
    params.input

    script:
    """
    mkdir ${model}
    Rscript ${baseDir}/bin/QQplot.R ${key}.txt ${model}/${key}
    """ 
}

process "GO_analysis" {

    label "low"
    label "finish"
    tag "${model}:${key}"
     
    input:
    tuple model, key, type,path("${model}/${key}.filtered_${params.output_FDR}_FDR.txt"},path(GOA),path(species)
    
    output:
    tuple model, key, type, path("GOA/BP_${model}/${key}.filtered_${params.output_FDR}/BP.txt"),path("GOA/MF_${model}/${key}.filtered_${params.output_FDR}/MF.txt"),path("GOA/CC_${model}/${key}.filtered_${params.output_FDR}/CC.txt")

    when:
    params.GOA

    script:
    """
    mkdir BP_${model}/${key}.filtered_${params.GO_filter}
    mkdir MF_${model}/${key}.filtered_${params.GO_filter}
    mkdir CC_${model}/${key}.filtered_${params.GO_filter}
    
    awk -F":" '\$1=\$1' ${model}/${key}.filtered_${params.GO_filter}_FDR.txt | awk -F"-" '\$1=\$1' | awk '{print \$1"\\t"\$2"\\t"\$3}' | sed '1d' > 2${model}/${key}.filtered_${params.GO_filter}_FDR.txt
    bedtools intersect -a ${GOA} -b 2${model}/${key}.filtered_${params.GO_filter}_FDR.txt | awk '\$3=="gene"' | awk -F";" '\$1=\$1' | awk '{gsub(/\ID=/,"",$9)}' | awk '{print \$9}' > 3${model}/${key}.filtered_${params.GO_filter}_FDR.txt
        
    bash ${baseDir}bin/GOA/GOtest.sh 3${model}/${key}.filtered_${params.GO_filter}_FDR.txt ${species} BP_${model}/${key}.filtered_${params.GO_filter}/overtree BP_${model}/${key}.filtered_${params.GO_filter}/over BP_${model}/${key}.filtered_${params.GO_filter}/undertree BP_${model}/${key}.filtered_${params.GO_filter}/under BP ${params.GO_filter}
    bash ${baseDir}bin/GOA/GOtest.sh 3${model}/${key}.filtered_${params.GO_filter}_FDR.txt ${species} MF_${model}/${key}.filtered_${params.GO_filter}/overtree MF_${model}/${key}.filtered_${params.GO_filter}/over MF_${model}/${key}.filtered_${params.GO_filter}/undertree MF_${model}/${key}.filtered_${params.GO_filter}/under MF ${params.GO_filter}
    bash ${baseDir}bin/GOA/GOtest.sh 3${model}/${key}.filtered_${params.GO_filter}_FDR.txt ${species} CC_${model}/${key}.filtered_${params.GO_filter}/overtree CC_${model}/${key}.filtered_${params.GO_filter}/over CC_${model}/${key}.filtered_${params.GO_filter}/undertree CC_${model}/${key}.filtered_${params.GO_filter}/under CC ${params.GO_filter}
    cat BP_${model}/${key}.filtered_${params.GO_filter}/*.txt > BP.txt
    cat MF_${model}/${key}.filtered_${params.GO_filter}/*.txt > MF.txt
    cat BCC_${model}/${key}.filtered_${params.GO_filter}/*.txt > CC.txt
    """
    
}
