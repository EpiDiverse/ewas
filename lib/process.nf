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


// GEM_Gmodel.out.map{ tuple( it[0] + "." + it[1], *it) }.groupTuple().map{ tuple("Gmodel", *it) }
// process to calculate FDR on combined files after splitting
process "calculate_FDR" {

    label "low"
    label "finish"
    tag "${model} - ${key}.txt"
     
    input:
    tuple model, key, contexts, types, path("txt"), path(results)
    
    output:
    tuple model, key, val("${types.unique().join("")}"), path("${model}/*.txt")
    tuple model, key, val("${types.unique().join("")}"), path(txt), path("${model}/${key}.filtered_${params.output_FDR}_FDR.txt")

    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel || params.GxE)

    script:
    """
    mkdir input ${model}
    tail -q -n+2 ${results} > input/${key}.txt

    if [[ \$(cat input/${key}.txt | wc -l) == 0 ]]; then
    echo "No findings with ${model == "Gmodel" ? "--Gmodel_pv ${params.Gmodel_pv}" : "--GxE_pv ${params.GxE_pv}"}" > ${model}/${key}.txt
    else
    Rscript ${baseDir}/bin/FDR.R input/${key}.txt ${model}/${key}.txt || exit \$?
    awk '{if(NR==1){print} else {if(\$6<=${params.output_FDR}){print}}}' ${model}/${key}.txt > ${model}/${key}.filtered_${params.output_FDR}_FDR.txt
    fi
    """ 
}



// GEM_Emodel.out[0]
// process to generate manhattan plots from Emodel
process "manhattan" {

    label "low"
    label "ignore"
    tag "${context} - ${type}"
     
    input:
    tuple context, type, path(txt)
    
    output:
    tuple type, path("*.png") optional true
    tuple type, path("*.zip") optional true

    when:
    params.input && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Emodel)

    script:
    """
    awk -F "\\t" 'BEGIN{print "SNP","CHR","BP","P"} NR!=1{split(\$1,cpg,":"); split(cpg[2],cpos,"-"); pos=(cpos[1]+cpos[2])/2;
    print \$1,cpg[1],pos,\$5}' ${txt} > manhattan.txt
    
    Rscript ${baseDir}/bin/manhattan.R manhattan.txt \$(basename ${txt} .txt) 0.00000001
    """ 
}



// calculate_FDR.out[0].filter{ it[0] == "Gmodel" }
// process to generate dotplots from Gmodel
process "dotPlot" {

    label "low"
    label "ignore"
    tag "${key}.txt"
     
    input:
    tuple model, key, type, path(result)
    
    output:
    tuple type, path("${model}/*.png") optional true

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
    tag "${key}.txt"
     
    input:
    tuple model, key, type, path("txt"), path(result)
    path snp
    path gxe
    
    output:
    tuple type, path("${model}/${key}/*.png")

    when:
    params.kplots > 0

    script:
    """
    mkdir ${model} ${model}/${key}
    head -qn 1 txt* | uniq > ${model}/${key}.txt
    tail -q -n+2 txt* >> ${model}/${key}.txt || exit \$?

    head -1 ${result} > ${model}/${key}/${result}
    tail -n+2 ${result} | sort -gk6 | head -${params.kplots} >> ${model}/${key}/${result} || exit \$?

    Rscript ${baseDir}/bin/Kplot.R ${model}/${key}/${result} ${model}/${key}.txt ${snp} ${gxe}
    """ 
}

