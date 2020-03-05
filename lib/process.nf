#!/usr/bin/env nextflow

// parse the samples file, separate the columns, and transpose
process "parsing" {

    label "low"
    label "finish"
    tag "${params.samples}"

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
    tuple model, key, contexts, types, path(txt), path(results)
    
    output:
    tuple model, key, val("${types.unique().join("")}"), path(txt), path("${model}/${key}.txt")
    tuple val("${contexts.unique().join("")}"), path("${model}/*.txt")

    when:
    params.input

    script:
    """
    mkdir input ${model}
    tail -q -n+2 ${results} > input/${key}.txt
    Rscript ${baseDir}/bin/FDR.R input/${key}.txt ${model}/${key}.txt
    awk '{if(NR==1){print} else {if(\$6<=${params.output_FDR}){print}}}' ${model}/${key}.txt > ${model}/${key}.filtered_${params.output_FDR}_FDR.txt
    """ 
}

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


// emodel_channel = bedtools_unionbedg.out.filter{it[1] == "bedGraph"}.mix(bedtools_intersect.out, average_over_regions.out)
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
    params.input
    
    script: 
    """
    awk -F "\\t" '{printf \"%s:%s-%s\",\$1,\$2,\$3; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' ${meth} > ${context}.txt
    Rscript ${baseDir}/bin/GEM_Emodel.R ${envs} ${covs} ${context}.txt ${params.Emodel_pv} ${context}.${type} > ${context}.${type}.log
    sort -V ${context}.${type}.txt |
    awk 'BEGIN{OFS="\\t"; print "cpg\\tbeta\\tstats\\tpvalue\\tFDR"} \$5<=${params.FDR}{print \$1,\$2,\$3,\$4,\$5}' > ${context}.${type}.filtered_${params.FDR}_FDR.txt
    """
}



/*

// something for SNPs
process "process_SNP_pipeline_SNP_file" {

    label "low"
    if(!params.DMP && !params.DMR){publishDir "${output_path}/meth_calls", mode:'copy'}
    if(params.DMP){publishDir "${output_path}/DMPs", mode: 'copy'}
    if(params.DMR){publishDir "${output_path}/DMRs", mode: 'copy'}
     
    input:
    file(snps) from snp1
    
    output:
    file("snps.txt") into (snp_methcalls, snp_DMPs, snp_DMRs)
    file("missing_stats.log")
    file("out.imiss")
    file("out.log")

    when:
    params.GEM_Gmodel && params.SNP

    script:
    """
    #!/bin/bash
    
    for i in *.vcf; do name=\$(basename "\$i" .vcf) ; sed -i "s|no_sample_specified|\$name|g" \$i ; bgzip *.vcf ; done
    for F in *.vcf.gz ; do   tabix -f -p vcf \${F}  ; done
    bcftools merge *vcf.gz -Oz -o merged.vcf.gz
    bcftools norm -Ov -m-snps merged.vcf.gz > splitted_merged.vcf.gz
    vcftools --gzvcf splitted_merged.vcf.gz --max-missing ${params.max_missing} --out missing_stats
    vcftools --gzvcf splitted_merged.vcf.gz --max-missing ${params.max_missing} --mac ${params.mac} --minQ ${params.minQ} --extract-FORMAT-info GT --out splitted_requiredFormat
    awk '{\$1=\$1"_"\$2; \$2=""; print \$0}' splitted_requiredFormat.GT.FORMAT | sed 's/^/SNP_/' | perl -pe 's!0/0!1!g; s!0/1!2!g; s!1/0!2!g; s!1/1!3!g; s!./.!0!g; s/CHROM_POS/ID/g' | awk -v OFS="\t" '\$1=\$1' |  awk 'NR == 1; NR > 1 {print \$0 | "sort -n"}' | uniq > snps.txt
    """  
} 

// something for SNPs
process "process_user_SNP_file" {

    label "low"
    
    if(!params.DMP && !params.DMR){publishDir "${output_path}/meth_calls", mode:'copy'}
    if(params.DMP){publishDir "${output_path}/DMPs", mode: 'copy'}
    if(params.DMR){publishDir "${output_path}/DMRs", mode: 'copy'}
     
    input:
    file snps from my_snp
    
    output:
    
    file("my_snps.txt") into (my_snp_methcalls, my_snp_DMPs, my_snp_DMRs)
    file("missing_stats.log")
    file("out.imiss")
    file("out.log")

    when:
    params.GEM_Gmodel && params.my_SNP

    script:
    """
    #!/bin/bash
    bcftools norm -Ov -m-snps merged.vcf.gz > splitted_merged.vcf.gz
    vcftools --gzvcf splitted_merged.vcf.gz --max-missing ${params.max_missing} --out missing_stats
    vcftools --gzvcf splitted_merged.vcf.gz --max-missing ${params.max_missing} --mac ${params.mac} --minQ ${params.minQ} --extract-FORMAT-info GT --out splitted_requiredFormat
    awk '{\$1=\$1"_"\$2; \$2=""; print \$0}' splitted_requiredFormat.GT.FORMAT | sed 's/^/SNP_/' | perl -pe 's!0/0!1!g; s!0/1!2!g; s!1/0!2!g; s!1/1!3!g; s!./.!0!g; s/CHROM_POS/ID/g' | awk -v OFS="\t" '\$1=\$1' |  awk 'NR == 1; NR > 1 {print \$0 | "sort -n"}' | uniq > my_snps.txt
    """  
} 

// eg. perform samtools index + make an additional file
process "process2" {

    label "low"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    input:
    tuple type, sample, path(bam)
    // eg. [process2, sample, /path/to/sorted.bam]

    output:
    tuple type, sample, path("${sample}.${type}.bam"), path("${sample}.${type}.bam.bai")
    // eg. [process2, sample, /path/to/process2.bam, /path/to/process2.bam.bai]
    path "${sample}.txt"

    // a directive when you only want the process to execute under certain conditions
    when:
    params.process2

    script:
    """
    cp sorted.bam ${sample}.${type}.bam
    samtools index ${sample}.${type}.bam
    touch ${sample}.txt
    """
}

// eg. perform samtools fastq
process "process3" {

    label "low"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    input:
    tuple type, sample, path(bam), path(bai)
    // eg. [process2, sample, /path/to/process2.bam, /path/to/process2.bam.bai]

    output:
    tuple sample, path("*.fastq.gz")
    // eg. [sample, /path/to/type.fastq.gz]
    // eg. [sample, [/path/to/type_1.fastq.gz, /path/to/type_2.fastq.gz]]

    script:
    if (params.PE)
        """
        samtools fastq -c 6 -1 ${sample}.${type}_1.fastq.gz -2 ${sample}.${type}_2.fastq.gz -0 /dev/null -s /dev/null -n ${bam}
        """
    else
        """
        samtools fastq -c 6 -0 /dev/null ${bam} > ${sample}.${type}.fastq.gz
        """
}


*/
=======

// calculate_FDR.out.filter{ it[0] == "GxEmodel" }
// process to generate top X plots from GxEmodel based on number provided by --kplots
process "topKplots" {

    label "low"
    label "finish"
    tag "${model} - ${key}.txt"
     
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
    head -1 txt > ${model}/${key}.txt
    tail -q -n+2 txt* >> ${model}/${key}.txt || exit \$?

    head -1 ${result} > ${model}/${key}/${result}
    tail -n+2 ${result} | sort -gk6 | head -${params.kplots} >> ${model}/${key}/${result} || exit \$?

    Rscript ${baseDir}/bin/kplot.R ${model}/${key}/${result} input/${key}.txt ${snp} ${gxe}
    """ 
}
