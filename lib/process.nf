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

    label "low"
    label "finish"
    tag "${model}:${key}"
     
    input:
    tuple model, key, contexts, types, path(results), path(logs)
    
    output:
    tuple model, key, val("${types.unique().join("")}"), path("${model}/*.txt")
    //tuple model, key, val("${types.unique().join("")}"), path("input/*.txt"), path("${model}/${key}.filtered_${params.output_FDR}_FDR.txt")

    when:
    params.input

    script:
    """
    mkdir input ${model}
    tail -q -n+2 ${results} > input/${key}.txt
    total=\$(cat ${logs} | grep "100.00%" | cut -d " " -f3 | tr -d "," | awk 'BEGIN{c=0} {c+=\$0} END{print c}')
    echo -e "${model == "Emodel" ? "cpg" : "cpg\\tsnp"}\\tbeta\\tstats\\tpvalue\\tFDR" |
    tee input/header.txt ${model}/${key}.txt ${model}/${key}.filtered_${params.output_FDR}_FDR.txt

    if [[ \$(head input/${key}.txt | wc -l) == 0 ]]; then
    echo "No findings with ${model == "Emodel" ? "--Emodel_pv ${params.Emodel_pv}" : model == "Gmodel" ? "--Gmodel_pv ${params.Gmodel_pv}" : "--GxE_pv ${params.GxE_pv}"}" > ${model}/${key}.txt
    else
    sort -grk${model == "Emodel" ? "4" : "5"} input/${key}.txt | cut -f${model == "Emodel" ? "2-" : "1-"} |
    awk -F "\\t" -v t="\$total" 'BEGIN{OFS="\\t";p=1;r=t} {fdr=(t/r)*${model == "Emodel" ? "\$4" : "\$5"};
    if(fdr>p){fdr=p}; if(fdr<=${params.output_FDR}){print \$0,fdr >> "${model}/${key}.filtered_${params.output_FDR}_FDR.txt"};
    print \$0,fdr; p=fdr;r--}' >> ${model}/${key}.txt || exit \$?
    fi
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
    awk -F "\\t" 'BEGIN{print "SNP","CHR","BP","P"} NR!=1{split(\$1,cpg,":"); split(cpg[2],cpos,"-"); pos=(cpos[1]+cpos[2])/2;
    print \$1,cpg[1],pos,\$5}' ${key}.filtered_${params.output_FDR}_FDR.txt > manhattan.txt
    Rscript ${baseDir}/bin/manhattan.R manhattan.txt ${key}.filtered_${params.output_FDR}_FDR 0.00000001 0.000001
    """ 
}



// calculate_FDR.out[0].filter{ it[0] == "Gmodel" }
// process to generate dotplots from Gmodel
process "dotPlot" {

    label "low"
    label "ignore"
    tag "${key}"
     
    input:
    tuple model, key, type, path(result)
    
    output:
    tuple type, path("${model}/*.png") optional true
    tuple type, path("${model}/*.zip") optional true

    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel)

    script:
    """
    mkdir ${model}
    awk -F "\\t" 'function abs(x){return ((x < 0.0) ? -x : x)}
    {if(NR!=1 && \$6<=${params.output_FDR}) {split(\$1,cpg,":"); split(\$2,snp,":"); split(cpg[2],cpos,"-"); split(snp[2],spos,"-");
    c=(cpos[1]+cpos[2])/2; s=(spos[1]+spos[2])/2;
    if(cpg[1]!=snp[1]){d="trans"} else {if(abs(c-s)>${params.distance}){d="trans"} else {d="cis"}};
    print cpg[1],c,snp[1],s,d}}' ${key}.filtered_${params.output_FDR}_FDR.txt > ${model}/${key}.txt
    
    Rscript ${baseDir}/bin/dotplot.R ${model}/${key}.txt ${model}/${key}.filtered_${params.output_FDR}_FDR 10
    """ 
}


// calculate_FDR.out[1].filter{ it[0] == "GxE" }
// process to generate top X plots from GxEmodel based on number provided by --kplots
process "topKplots" {

    label "low"
    label "ignore"
    tag "${key}"
     
    input:
    tuple key, type, path(result), path("txt")
    // eg. [CHG.region, region, [path/to/CHG.region.txt, /path/to/filtered.txt], [path/to/CHG.txt, path/to/CHG.txt, ...]]
    path snp
    path gxe
    
    output:
    tuple type, path("GxE/${key}/*.png") optional true

    when:
    params.kplots > 0

    script:
    """
    mkdir GxE GxE/${key}
    head -qn 1 txt* | uniq > GxE/${key}.txt
    tail -qn+2 txt* >> GxE/${key}.txt

    awk 'NR==1{print;next}{print | "sort -gk6 | head -${params.kplots}"}' ${key}.filtered_${params.output_FDR}_FDR.txt \\
    > GxE/${key}/${key}.filtered_${params.output_FDR}_FDR.txt || exit \$?
    Rscript ${baseDir}/bin/Kplot.R GxE/${key}/${key}.filtered_${params.output_FDR}_FDR.txt GxE/${key}.txt ${snp} ${gxe}
    """ 
}

