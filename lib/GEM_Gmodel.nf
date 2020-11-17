#!/usr/bin/env nextflow

// preparing single-sample vcf files for merging 
process "tabix" {

    label "low"
    label "finish"
    tag "$sample"
     
    input:
    tuple val(sample), path(snp)
    
    output:
    path "output/*.vcf.gz"
    path "output/*.tbi"
    
    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE && !params.GWAS) || params.Gmodel || params.GxE || params.GWAS)

    script:
    """
    mkdir output
    sed "s|no_sample_specified|${sample}|g" <(bcftools view ${snp}) > output/${sample}.vcf
    bgzip output/${sample}.vcf
    tabix -f -p vcf output/${sample}.vcf.gz
    """
}


//nextflow.processor.TaskPath
//nextflow.util.BlankSeparatedList

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
     params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE && !params.GWAS) || params.Gmodel || params.GxE || params.GWAS)

    script:
    """
    mkdir input
    ${snps.getClass() == nextflow.util.BlankSeparatedList && snps.size() > 1 ? "bcftools merge ${snps} -Oz -o input/merged.vcf.gz || exit \$?" : ""}
    bcftools norm -Ov -m-snps ${snps.getClass() == nextflow.util.BlankSeparatedList && snps.size() > 1 ? "input/merged.vcf.gz" : "${snps}"} > input/norm.vcf.gz || exit \$?

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
    path "out.recode.vcf.gz"
    //file("out.imiss")
    //file("out.log")
    //path "out.log"

    when:
     params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE && !params.GWAS) || params.Gmodel || params.GxE || params.GWAS)

    script:
    """
    vcftools --gzvcf ${snp} --max-missing ${params.max_missing} --recode --stdout | bgzip  > out.recode.vcf.gz
    """ 
} 

process "BEAGLE_SNP_Imputation" {
    
    label "low"
    label "finish"
    
    //beforeScript "set +u; source activate ewas"
    //afterScript "set +u; source deactivate"
    
    input:
    path "out.recode.vcf.gz"
    
    output:
    path "snps_imputed.gt.vcf.gz"

    when:
     params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE && !params.GWAS) || params.Gmodel || params.GxE || params.GWAS)
    
    script:
    """
    beagle gt=out.recode.vcf.gz iterations=${params.iters} phase-states=${params.phase_states} imp-states=${params.imp_states} ne=${params.ne} nthreads=${params.nthreads_SNP} out=snps_imputed.gt
   
    """ 
} 

// bcftools.out
// extract required format
process "vcftools_extract" {

    label "low"
    label "finish"
    
    publishDir "${params.output}/input", pattern: "snps.txt" , mode: 'copy', \
            enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel || params.GxE) ? true : false
    
    input:
    path ("snps_imputed.gt.vcf.gz")
    
    output:
    path "snps.txt"
    //file("out.imiss")
    //file("out.log")

    when:
     params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE && !params.GWAS) || params.Gmodel || params.GxE || params.GWAS)

    script:
    """   
    vcftools --gzvcf snps_imputed.gt.vcf.gz \\
    --max-missing 1 \\
    --mac ${params.mac} \\
    --minQ ${params.minQ} \\
    --012 \\
    --out GT || exit \$?

    paste <(cat <(echo -e "CHROM\\tPOS") GT.012.pos) <(paste GT.012.indv <(cut -f2- GT.012) | datamash transpose) |
    awk '{printf "%s:%s-%s",\$1,\$2-1,\$2; for(i=3; i<=NF; i++) {printf "\\t%s",(NR==1?\$i:\$i+1)}; print null}' > snps.txt
    """ 
}


//split_scaffolds.out.transpose()
// RUN GEM Gmodel
process "GEM_Gmodel" {
    
    label "low"
    label "finish"
    tag "${context}.${type} - ${meth.baseName}"
    
    //publishDir "${params.output}/positions/Gmodel", pattern: "*.txt" , mode: 'copy', \
    //        enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel) ? true : false
    //publishDir "${params.output}/regions/Gmodel", pattern: "*.txt" , mode: 'copy', \
    //        enabled: (params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel)) && params.DMRs ? true : false
    
    input:
    tuple val(context), val(type), path(meth)
    path snps
    path covs
    
    output:
    //tuple context, type, path("output/*.txt"), path("output/*.log")
    path "output/${context}.${type}.txt"
    path "output/${context}.${type}.log"
   
    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE && !params.GWAS) || params.Gmodel)
    
    script: 
    """
    mkdir output
    awk -F "\\t" '{printf \"%s:%s-%s\",\$1,\$2,\$3; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' ${meth} > \$(basename ${meth} .bed).txt
    Rscript ${baseDir}/bin/GEM_Gmodel.R ${baseDir}/bin ${snps} ${covs} \$(basename ${meth} .bed).txt ${params.Gmodel_pv} output/temp > output/${context}.${type}.log || exit \$?
    tail -n+2 output/temp.txt > output/${context}.${type}.txt
    #Rscript ${baseDir}/bin/GEM_Gmodel.R ${baseDir}/bin ${snps} ${covs} \$(basename ${meth} .bed).txt ${params.Gmodel_pv} output/\$(basename ${meth} .bed) > output/\$(basename ${meth} .bed).log
    """
}


// RUN GEM Gmodel
process "GEM_GxEmodel" {
    
    label "low"
    label "finish"
    tag "${context}.${type} - ${meth.baseName}"
    
    //publishDir "${params.output}/positions/GxE", pattern: "*.txt" , mode: 'copy', \
    //        enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) ? true : false
    //publishDir "${params.output}/regions/GxE", pattern: "*.txt" , mode: 'copy', \
    //        enabled: (params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE)) && params.DMRs ? true : false
    
    input:
    tuple val(context), val(type), path(meth)
    path snps
    path gxe
    
    output:
    //tuple context, type, path("output/*.txt"), path("output/*.log")
    //tuple context, type, path("*.txt")
    path "output/${context}.${type}.txt"
    path "output/${context}.${type}.log"
    path "${context}.${type}.txt"
    tuple val(context), val(type), path("header.txt")

    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE && !params.GWAS) || params.GxE)
    
    script: 
    """
    mkdir output
    
    #awk -F "\\t" '{if(NR==1){
    #printf "%s:%s-%s",\$1,\$2,\$3 > "header.txt"; printf "%s:%s-%s",\$1,\$2,\$3 > "meth.txt";
    #for(i=4; i<=NF; i++) {printf "\\t%s",\$i > "header.txt"; printf "\\t%s",\$i > "meth.txt"};
    #print null > "header.txt"; print null > "meth.txt"}else{
    #printf \"%s:%s-%s\",\$1,\$2,\$3; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}}' ${meth} |
    #tee -a meth.txt > ${context}.${type}.txt

    #awk -F "\\t" '{tee=(NR==1?"tee output/header.txt":"tee ${context}.${type}.txt");
    #printf "%s:%s-%s",\$1,\$2,\$3 | tee; for(i=4; i<=NF; i++) {printf "\\t%s",\$i | tee};
    #print null | tee}' ${meth} > meth.txt

    awk -F "\\t" '{printf \"%s:%s-%s\",\$1,\$2,\$3; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' ${meth} > meth.txt
    head -1 meth.txt > header.txt && tail -n+2 meth.txt > ${context}.${type}.txt
    #awk -F "\\t" '{printf \"%s:%s-%s\",\$1,\$2,\$3; for(i=4; i<=NF; i++) {printf \"\\t%s\",\$i}; print null}' ${meth} > \$(basename ${meth} .bed).txt
    Rscript ${baseDir}/bin/GEM_GxE.R ${baseDir}/bin ${snps} ${gxe} meth.txt ${params.GxE_pv} output/meth > output/${context}.${type}.log || exit \$?
    tail -n+2 output/meth.txt > output/${context}.${type}.txt
    #Rscript ${baseDir}/bin/GEM_GxE.R ${baseDir}/bin ${snps} ${gxe} \$(basename ${meth} .bed).txt ${params.GxE_pv} output/\$(basename ${meth} .bed) > output/\$(basename ${meth} .bed).log
    """
}


process "GEM_GWAS" {
    
    label "low"
    label "finish"
    tag "${context}.${type} - ${meth.baseName}"


    //enabled: (params.SNPs && ((!params.Emodel && !params.GWAS && !params.GxE) || params.GWAS)) && params.DMRs ? true : false
    
    input:
    path envs
    path snps
    path covs
    
    output:
    path "output/${context}.${type}.gz"
    path "output/${context}.${type}.log"
   
    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE && !params.GWAS) || params.GWAS)
    //params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel)
    
    script: 
    """
    mkdir output
    Rscript ${baseDir}/bin/GEM_GWASmodel.R ${baseDir}/bin ${snps} ${covs} ${envs} ${params.GWAS_pv} output/temp > output/${context}.${type}.log || exit \$?
    tail -n+2 output/temp.txt | awk 'BEGIN{OFS=\"\\t\"} {printf \"%s\\t%s\\t%s\\t%s\\t%s\\n\", \$2,\$1,\$3,\$4,\$5}' | gzip > output/${context}.${type}.gz  && rm output/temp.txt   
    """
}



// calculate_FDR.out[0].filter{ it[0] == "Gmodel" }
// process to generate dotplots from Gmodel
process "dotPlot" {

    label "low"
    label "ignore"
    tag "${key}"
    
    publishDir "${params.output}/positions", pattern: "${model}/*.bedGraph.filtered_${params.output_FDR}_FDR.png" , mode: 'copy', enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel) ? true : false
    publishDir "${params.output}/positions", pattern: "${model}/*.bedGraph.filtered_${params.output_FDR}_FDR.zip" , mode: 'copy', enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel) ? true : false    
    publishDir "${params.output}/positions", pattern: "${model}/*.DMRs.filtered_${params.output_FDR}_FDR.png" , mode: 'copy', enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel) ? true : false
    publishDir "${params.output}/positions", pattern: "${model}/*.DMRs.filtered_${params.output_FDR}_FDR.zip" , mode: 'copy', enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel) ? true : false    
    publishDir "${params.output}/positions", pattern: "${model}/*.DMPs.filtered_${params.output_FDR}_FDR.png" , mode: 'copy', enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel) ? true : false
    publishDir "${params.output}/positions", pattern: "${model}/*.DMPs.filtered_${params.output_FDR}_FDR.zip" , mode: 'copy', enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel) ? true : false    
    publishDir "${params.output}/regions", pattern: "${model}/*.region.filtered_${params.output_FDR}_FDR.png" , mode: 'copy', enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel) ? true : false
    publishDir "${params.output}/regions", pattern: "${model}/*.region.filtered_${params.output_FDR}_FDR.zip" , mode: 'copy', enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.Gmodel) ? true : false    
   
    
    input:
    tuple val(model), val(key), val(type), path(result)
    
    output:
    tuple val(type), path("${model}/*.png") optional true
    tuple val(type), path("${model}/*.zip") optional true

    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE && !params.GWAS) || params.Gmodel)

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
    
    publishDir "${params.output}/regions", pattern: "${model}/CpG.region/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
    
    publishDir "${params.output}/regions", pattern: "${model}/CHG.region/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
    
    publishDir "${params.output}/regions", pattern: "${model}/CHH.region/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
    
    publishDir "${params.output}/positions", pattern: "${model}/CpG.bedGraph/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
    
    publishDir "${params.output}/positions", pattern: "${model}/CHG.bedGraph/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
    
    publishDir "${params.output}/positions", pattern: "${model}/CHH.bedGraph/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
  
    publishDir "${params.output}/positions", pattern: "${model}/CpG.DMRs/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
    
    publishDir "${params.output}/positions", pattern: "${model}/CHG.DMRs/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
    
    publishDir "${params.output}/positions", pattern: "${model}/CHH.DMRs/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
 
 
    publishDir "${params.output}/positions", pattern: "${model}/CpG.DMPs/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
    
    publishDir "${params.output}/positions", pattern: "${model}/CHG.DMPs/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
    
    publishDir "${params.output}/positions", pattern: "${model}/CHH.DMPs/*.png" , mode: 'copy', \
    enabled: params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0 ? true : false
    
    
    input:
    //tuple key, type, path(result), path(scaffolds)
    // eg. [CHG.region, region, [path/to/CHG.region.txt, /path/to/filtered.txt], path/to/CHG.region.txt]
    tuple val(key), val(type), path(results), path("meth.txt"), path(header)
    path snp
    path gxe
    
    output:
    tuple val(type), path("GxE/${key}/*.png") optional true

    when:
    params.SNPs && ((!params.Emodel && !params.Gmodel && !params.GxE) || params.GxE) && params.kplots > 0

    script:
    """
    mkdir GxE GxE/${key}
    #head -n 1 \${scaffolds} > GxE/${key}.txt
    #tail -n+2 \${scaffolds} >> GxE/${key}.txt

    awk 'NR==1{print;next}{print | "sort -gk6 | head -${params.kplots}"}' ${key}.filtered_${params.output_FDR}_FDR.txt \\
    > GxE/${key}/${key}.filtered_${params.output_FDR}_FDR.txt || exit \$?
    Rscript ${baseDir}/bin/Kplot.R GxE/${key}/${key}.filtered_${params.output_FDR}_FDR.txt <(cat ${header} meth.txt) ${snp} ${gxe}
    """ 
}
