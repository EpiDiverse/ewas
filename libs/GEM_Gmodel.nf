#!/usr/bin/env nextflow

// preparing single-sample vcf files for merging 
process "tabix" {

    label "low"
    label "finish"
    tag "$sample"
     
    input:
    tuple sample, path(snp)
    
    output:
    path "output/*.vcf.gz"
    path "output/*.tbi"
    
    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel)

    script:
    """
    mkdir output
    sed "s|no_sample_specified|${sample}|g" <(bcftools view ${snp}) > output/${sample}.vcf
    bgzip output/${sample}.vcf
    tabix -f -p vcf output/${sample}.vcf.gz
    """
}


// tabix.out[0].collect()
// tabix.out[1].collect()
// merge single-sample vcf files
process "bcftools" {

    label "low"
    label "finish"
     
    input:
    path samples
    path snps
    path tbis
    
    output:
    path "input/filtered.vcf.gz"

    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel)

    script:
    """
    mkdir input
    ${snps.size() > 1 ? "bcftools merge ${snps} -Oz -o input/merged.vcf.gz || exit \$?" : ""}
    bcftools norm -Ov -m-snps ${snps.size() > 1 ? "input/merged.vcf.gz" : "${snps}"} > input/norm.vcf.gz || exit \$?

    bcftools view -S <(cut -f1 ${samples}) input/norm.vcf.gz > input/filtered.vcf.gz || exit \$?
    bcftools query -l input/filtered.vcf.gz > input/samples.txt || exit \$?

    total=\$(cat ${samples} | wc -l)
    match=\$(grep -f <(cut -f1 ${samples}) input/samples.txt | wc -l)

    if [[ \$match != \$total ]];
    then echo "ERROR: multi-sample vcf is missing sample IDs from: ${samples}"; exit 1;
    fi
    """
}


// bcftools.out
// get missing information
process "vcftools_missing" {

    label "low"
    label "finish"
     
    input:
    path snp
    
    output:
    path "missing_stats.log"
    //file("out.imiss")
    //file("out.log")
    //path "out.log"

    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel)

    script:
    """
    vcftools --gzvcf ${snp} --max-missing ${params.max_missing} > missing_stats.log
    """ 
} 


// bcftools.out
// extract required format
process "vcftools_extract" {

    label "low"
    label "finish"
     
    input:
    path snp
    
    output:
    path "snps.txt"
    //file("out.imiss")
    //file("out.log")

    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel)

    script:
    """   
    vcftools --gzvcf ${snp} \\
    --max-missing ${params.max_missing} \\
    --mac ${params.mac} \\
    --minQ ${params.minQ} \\
    --extract-FORMAT-info GT \\
    --out splitted_requiredFormat || exit \$?

    awk '{\$1=\$1"_"\$2; \$2=""; print \$0}' splitted_requiredFormat.GT.FORMAT |
    sed 's/^/SNP_/' |
    perl -pe 's!0/0!1!g; s!0/1!2!g; s!1/0!2!g; s!1/1!3!g; s!./.!0!g; s/CHROM_POS/ID/g' |
    awk -v OFS="\t" '\$1=\$1' |
    awk 'NR == 1; NR > 1 {print \$0 | "sort -n"}' |
    uniq > snps.txt
    """ 
}

// meth_channel = bedtools_unionbedg.out.filter{it[1] == "bedGraph"}.mix(bedtools_intersect.out, average_over_regions.out)
// eg. [CpG, bedGraph, methylation.txt]
// eg. [CpG, DMRs, CpG.bed]
// eg. [CpG, regions, CpG.avg]

// RUN GEM Gmodel
process "GEM_Gmodel" {
    
    label "low"
    label "finish"
    tag "$type - $context"

    input:
    tuple context, type, path(meth)
    path snps
    path covs
    
    output:
    tuple type, path("${context}.${type}.filtered_${params.FDR}_FDR.txt")
    tuple type, path("${context}.${type}.txt")
    tuple type, path("${context}.${type}.jpg")
    tuple type, path("${context}.${type}.log")
   
    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel)
    
    script: 
    """
    awk -F "\\t" '{printf \"%s:%s-%s\",\$1,\$2,\$3; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' ${meth} > ${context}.txt
    Rscript ${baseDir}/bin/GEM_Emodel.R ${snps} ${covs} ${context}.txt ${params.Gmodel_pv} ${context}.${type} > ${context}.${type}.log
    sort -V ${context}.${type}.txt |
    awk 'BEGIN{OFS="\\t"; print "cpg\\tsnp\\tbeta\\tstats\\tpvalue\\tFDR"} \$6<=${params.FDR}{print \$1,\$2,\$3,\$4,\$5,\$6}' > ${context}.${type}.filtered_${params.FDR}_FDR.txt
    """
}