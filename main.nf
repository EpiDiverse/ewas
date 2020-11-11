#!/usr/bin/env nextflow

// DSL2 BRANCH
nextflow.enable.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\

         =================================================
          E P I D I V E R S E - E W A S   P I P E L I N E
         =================================================
         ~ version ${workflow.manifest.version}

         Usage: 
              nextflow run epidiverse/ewas [OPTIONS]...

         Options: GENERAL
         
              --input [path/to/input/dir]      [REQUIRED] Specify input path for the directory containing outputs from the WGBS pipeline.
                                               The pipeline searches for bedGraph files in '*/bedGraph/{sample_name}_{context}.bedGraph'
                                               format, where sample names must correspond to the samplesheet and context can be either
                                               "CpG", "CHG", or "CHH".

              --samples [path/to/samples.tsv]  [REQUIRED] Specify the path to the samplesheet file containing information regarding
                                               sample names and corresponding environment and covariate values. The file must contain
                                               at least three tab-separated columns: 1) sample names, 2) environment value, 3) covariate
                                               values, with further columns optional for additional covariates.

              --DMPs [path/to/DMPs/dir]        Specify path to the DMR pipeline output directory to run EWAS analyses in addition with
                                               methylated positions filtered by significant DMPs. The pipeline searches for bed files in
                                               '*/{context}/metilene/*/*.bed' format where context can be either "CpG", "CHG", or "CHH".
         
              --DMRs [path/to/DMRs/dir]        Specify path to the DMR pipeline output directory to run EWAS analyses in addition with
                                               methylated positions filtered by significant DMRs. In addition, the pipeline will call
                                               the union of all significant regions and attempt to run EWAS with whole regions as markers.
                                               The pipeline searches for bed files in '*/{context}/metilene/*/*.bed' format where context
                                               can be either "CpG", "CHG", or "CHH".

              --SNPs [path/to/vcf/dir]         Specify path to the SNP pipeline output directory to enable EWAS analyses Gmodel and GxEmodel
                                               which attempt to create a genome-wide methQTL map. ONLY SUITABLE FOR DIPLOID ORGANISMS.
                                               The pipeline searches for VCF files in '*/vcf/{sample_name}.{extension}' where sample names
                                               must correspond to the samplesheet and the extension can be any standard vcf extension
                                               readable by 'bcftools' and defined with --extension parameter. Alternatively, the path to a
                                               single multi-sample VCF file can be provided.

              --extension <STR>                Specify the extension to use when searching for VCF files [default: vcf.gz]

              --output <STR>                   A string that can be given to name the output directory [default: ewas]


         Options: MODEL DECISION
              --Emodel                         Run analysis with "E model". Disables other models unless they are also specified. If no
                                               individual model is specified then all that are possible with the provided inputs will run
                                               in parallel [default: off]

              --Gmodel                         Run analysis with "G model". Disables other models unless they are also specified. If no
                                               individual model is specified then all that are possible with the provided inputs will run
                                               in parallel [default: off]

              --GxE                            Run analysis with "GxE model". Disables other models unless they are also specified. If no
                                               individual model is specified then all that are possible with the provided inputs will run
                                               in parallel [default: off]
            
              --noCpG                          Disables EWAS analysis in CpG context. Note: at least one methylation context is required
                                               for analysis. [default: off]

              --noCHG                          Disables EWAS analysis in CHG context. Note: at least one methylation context is required
                                               for analysis. [default: off]

              --noCHH                          Disables EWAS analysis in CHH context. Note: at least one methylation context is required
                                               for analysis. [default: off]

              --all                            If --DMPs and/or --DMRs are provided to the pipeline then raw un-intersected bedGraphs are
                                               not carried forward for analysis. Enable this parameter to process them alongside the 
                                               DMP / DMR intersections in parallel [default: off]

              --input [path/to/input/dir]       [REQUIRED] Specify input path for the directory containing outputs from the WGBS pipeline.
                                                The pipeline searches for bedGraph files in '*/bedGraph/{sample_name}_{context}.bedGraph'
                                                format, where sample names must correspond to the samplesheet and context can be either
                                                "CpG", "CHG", or "CHH".

              --samples [path/to/samples.tsv]   [REQUIRED] Specify the path to the samplesheet file containing information regarding
                                                sample names and corresponding environment and covariate values. The file must contain
                                                at least three tab-separated columns: 1) sample names, 2) environment value, 3) covariate
                                                values, with further columns optional for additional covariates.

              --DMPs [path/to/DMPs/dir]         Specify path to the DMR pipeline output directory to run EWAS analyses in addition with
                                                methylated positions filtered by significant DMPs. The pipeline searches for bed files in
                                                '*/{context}/metilene/*/*.bed' format where context can be either "CpG", "CHG", or "CHH".
         
              --DMRs [path/to/DMRs/dir]         Specify path to the DMR pipeline output directory to run EWAS analyses in addition with
                                                methylated positions filtered by significant DMRs. In addition, the pipeline will call
                                                the union of all significant regions and attempt to run EWAS with whole regions as markers.
                                                The pipeline searches for bed files in '*/{context}/metilene/*/*.bed' format where context
                                                can be either "CpG", "CHG", or "CHH".

              --SNPs [path/to/vcf/dir]          Specify path to the SNP pipeline output directory to enable EWAS analyses Gmodel and GxEmodel
                                                which attempt to create a genome-wide methQTL map. ONLY SUITABLE FOR DIPLOID ORGANISMS.
                                                The pipeline searches for VCF files in '*/vcf/{sample_name}.{extension}' where sample names
                                                must correspond to the samplesheet and the extension can be any standard vcf extension
                                                readable by 'bcftools' and defined with --extension parameter. Alternatively, the path to a
                                                single multi-sample VCF file can be provided.

              --extension <STR>                 Specify the extension to use when searching for VCF files [default: vcf.gz]

              --output <STR>                    A string that can be given to name the output directory [default: ewas]


         Options: MODEL DECISION
              --Emodel                          Run analysis with "E model". Disables other models unless they are also specified. If no
                                                individual model is specified then all that are possible with the provided inputs will run
                                                in parallel [default: off]

              --Gmodel                          Run analysis with "G model". Disables other models unless they are also specified. If no
                                                individual model is specified then all that are possible with the provided inputs will run
                                                in parallel [default: off]

              --GxE                             Run analysis with "GxE model". Disables other models unless they are also specified. If no
                                                individual model is specified then all that are possible with the provided inputs will run
                                                in parallel [default: off]
            
              --noCpG                           Disables EWAS analysis in CpG context. Note: at least one methylation context is required
                                                for analysis. [default: off]

              --noCHG                           Disables EWAS analysis in CHG context. Note: at least one methylation context is required
                                                for analysis. [default: off]

              --noCHH                           Disables EWAS analysis in CHH context. Note: at least one methylation context is required
                                                for analysis. [default: off]

              --all                             If --DMPs and/or --DMRs are provided to the pipeline then raw un-intersected bedGraphs are
                                                not carried forward for analysis. Enable this parameter to process them alongside the 
                                                DMP / DMR intersections in parallel [default: off]


         Options: INPUT FILTERING
              --coverage <INT>                 Specify the minimum coverage threshold to filter individual methylated positions from the
                                               --input directory before running analyses [default: 0] 
            
              --filter_FDR <FLOAT>             Specify the minimum FDR significance threshold to include DMPs and/or DMRs from the respective
                                               --DMPs and --DMRs directories [default: 0.05]          

              --filter_NA <FLOAT>              Specify the maximum proportion of samples that can contain a missing value before a methylated

                                               position is removed from the analysis Please be careful about the usage, filtering out more than 
                                               half of the data is not recommended.[default: 0] . 
                                               

                                               position is removed from the analysis [default: 0] 

              --filter_SD <FLOAT>              Specify the maximum standard deviation in methylation between samples to filter individual
                                               positions based on the degree of difference [default: 0] 

              --proportion <FLOAT>             Minimum proportion of samples that must share a DMP and/or DMR for it to be considered in the
                                               analysis [default: 0.2]

              --merge                          When running EWAS using the union set of DMRs as markers, specify to merge adjacent sub-regions
                                               into larger regions prior to methylation averaging and subsequent analysis [default: off]


         Options: SNP FILTERING
              --mac <INT>                      Minor allele count [default: 3]



              --minQ <INT>                     Minimum quality score [default: 30]


         Options: OUTPUT FILTERING       

              
              --minQ <INT>                     Minimum quality score [default: 30]


         Options: Gene Ontology(GO) Analysis
              --goa                            Path to .gff file with standard format that has name of chromosome/scaffold (the first column),
                                               name of program generated this feature (the second column), feature type name, e.g. Gene, Variation,
                                               Similarity (the third column), start (fourth column), end positions (fifth column), score, a floating
                                               point (the sixth column), strand (the seventh column), frame (the eighth column) and attribute (the 
                                               last column).
              
              --GO_filter <FLOAT>              Specify the minimum p-value for filtering the output [default: 0.05]
             
              --species                        Specify the name of the species to run GOA. Possible options: "fragaria_vesca", "picea_abies", 
                                               "populus_nigra" and "thlaspi_arvense".
                                               

         Options: OUTPUT FILTERING   

              --output_FDR <FLOAT>            Specify the maximum FDR threshold for filtering EWAS post-analysis [default: 0.05]
              
              --Emodel_pv <FLOAT>             Set the p-value to run "E model". Note: this filter is applied prior to FDR calculation


              --Gmodel_pv <FLOAT>             Set the p-value to run "G model". Note: this filter is applied prior to FDR calculation
                                              and should be used cautiously [default: 1]

              --GxE_pv <FLOAT>                Set the p-value to run "GxE model". Note: this filter is applied prior to FDR calculation
                                              and should be used cautiously [default: 1]


                                              
              --Gmodel_pv <FLOAT>             Set the p-value to run "G model". Note: this filter is applied prior to FDR calculation
                                              and should be used cautiously [default: 1]
                                              
              --GxE_pv <FLOAT>                Set the p-value to run "GxE model". Note: this filter is applied prior to FDR calculation
                                              and should be used cautiously [default: 1]


         Options: VISUALISATION
              --kplots <INT>                  Specify the number of plots to generate for the top k significant results in "GxE model"  
                                              [default: 10]

              --distance <INT>                Specify the distance threshold to define cis and trans methQTLs in the dotplot generated
                                              for "G model" output [default: 5000]

         Options: ADDITIONAL PARAMS
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    


         Example: 
              nextflow run epidiverse/ewas \
              --input /path/to/input/dir \
              --samples /path/to/samples.tsv
              --output ewas

    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// PRINT VERSION AND EXIT
if(params.version){
    println """\
         =================================================
          E P I D I V E R S E - E W A S   P I P E L I N E
         =================================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// VALIDATE PARAMETERS
ParameterChecks.checkParams(params)

// DEFINE PATHS
CpG_path = "${params.input}/*/bedGraph/*_CpG.bedGraph"
CHG_path = "${params.input}/*/bedGraph/*_CHG.bedGraph"
CHH_path = "${params.input}/*/bedGraph/*_CHH.bedGraph"
//CpG_path = "${params.input}/CpG/*.bedGraph"
//CHG_path = "${params.input}/CHG/*.bedGraph"
//CHH_path = "${params.input}/CHH/*.bedGraph"
  
CpG_path_DMPs = "${params.DMPs}/CpG/metilene/*" 
CHG_path_DMPs = "${params.DMPs}/CHG/metilene/*"
CHH_path_DMPs = "${params.DMPs}/CHH/metilene/*"
//CpG_path_DMPs = "${params.DMPs}/CpG/*.bed" 
//CHG_path_DMPs = "${params.DMPs}/CHG/*.bed"
//CHH_path_DMPs = "${params.DMPs}/CHH/*.bed"

CpG_path_DMRs = "${params.DMRs}/CpG/metilene/*" 
CHG_path_DMRs = "${params.DMRs}/CHG/metilene/*"
CHH_path_DMRs = "${params.DMRs}/CHH/metilene/*"
//CpG_path_DMRs = "${params.DMRs}/CpG/*.bed" 
//CHG_path_DMRs = "${params.DMRs}/CHG/*.bed"
//CHH_path_DMRs = "${params.DMRs}/CHH/*.bed"

// PARAMETER CHECKS
//if( !params.input ){error "ERROR: Missing required parameter --input"}
//if( params.noCpG && params.noCHG && params.noCHH ){error "ERROR: please specify at least one methylation context for analysis"}
if( !params.Emodel && !params.Gmodel && !params.GxE ){
    Emodel = true
    Gmodel = true
    GxE = true
} else {
    Emodel = params.Emodel
    Gmodel = params.Gmodel
    GxE = params.GxE
}


// PRINT STANDARD LOGGING INFO
log.info ""
log.info "         =================================================="
log.info "          E P I D I V E R S E - E W A S    P I P E L I N E"
if(params.debug){
log.info "         (debug mode enabled)"
log.info "         ==================================================" }
else {
log.info "         ==================================================" }
log.info "         ~ version ${workflow.manifest.version}"
log.info ""
log.info "         samplesheet                     : ${params.samples}"
log.info "         input dir                       : ${params.input}"
log.info "         DMPs dir                        : ${params.DMPs ? "$params.DMPs" : "-"}"
log.info "         DMRs dir                        : ${params.DMRs ? "$params.DMRs" : "-"}"
log.info "         VCF(s)                          : ${params.SNPs ? "$params.SNPs" : "-"}"
log.info "         output dir                      : ${params.output}"
log.info ""
log.info "         Analysis Configuration"
log.info "         =================================================="
log.info "         GEM model(s)                    : ${Emodel ? "Emodel " : ""}${Gmodel ? "Gmodel " : ""}${GxE ? "GxE" : ""}"
log.info "         Methylation context(s)          : ${params.noCpG ? "" : "CpG "}${params.noCHH ? "" : "CHH "}${params.noCHG ? "" : "CHG"}"
if(Emodel && params.Emodel_pv.toInteger() < 1){
log.info "         Emodel p-value                  : ${params.Emodel_pv}" }
if(params.SNPs && Gmodel && params.Gmodel_pv.toInteger() < 1){
log.info "         Gmodel p-value                  : ${params.Gmodel_pv}" }
if(params.SNPs && GxE && params.GxE_pv.toInteger() < 1){
log.info "         GxE p-value                     : ${params.GxE_pv}" }
log.info ""
log.info "         Input Filtering"
log.info "         =================================================="
log.info "         coverage                        : ${params.coverage}"
log.info "         NA filtering                    : ${params.filter_NA}"
log.info "         SD filtering                    : ${params.filter_SD}"
if(params.DMPs || params.DMRs){
log.info "         input FDR                       : ${params.filter_FDR}"
log.info "         overlap proportion              : ${params.proportion}" }
if(params.DMRs){
log.info "         merge regions?                  : ${params.merge ? "True" : "False"}" }
log.info ""
if(params.SNPs){
log.info "         =================================================="
log.info "         SNP Filtering"
log.info "         =================================================="
log.info "         maximum missing                 : ${params.max_missing}"
log.info "         minor allele count              : ${params.mac}"
log.info "         minimum quality score           : ${params.minQ}"
log.info "" }

if(params.burnin || params.iterations || params.phase_states || params.imp_states || params.ne || nthreads_SNP){
log.info "         =================================================="
log.info "         SNP Imputation with BEAGLE"
log.info "         =================================================="
log.info "         burnin iterations                                      : ${params.burnin}"
log.info "         iterations for genotype phase estimation               : ${params.iters}"
log.info "         number of model states for genotype estimation         : ${params.phase_states}"
log.info "         number of model states for ungenotype estimation       : ${params.imp_states}"
log.info "         effective population size                              : ${params.ne}"
log.info "         number of threads of execution                         : ${params.nthreads_SNP}"
log.info "" }

/*
if(params.goa && params.species || params.GO_filter){
log.info "         =================================================="
log.info "         Gene Ontology Analysis"
log.info "         =================================================="
log.info "         path to .gff file               : ${params.goa}"
log.info "         name of species                 : ${params.species}"
log.info "         output GO p-value               : ${params.GO_filter}"
log.info "" }
*/
//if(params.burnin || params.iterations || params.phase_states || params.imp_states || params.ne || nthreads_SNP){
//log.info "         =================================================="
//log.info "         SNP Imputation with BEAGLE"
//log.info "         =================================================="
//log.info "         burnin iterations                                      : ${params.burnin}"
//log.info "         iterations for genotype phase estimation               : ${params.iters}"
//log.info "         number of model states for genotype estimation         : ${params.phase_states}"
//log.info "         number of model states for ungenotype estimation       : ${params.imp_states}"
//log.info "         effective population size                              : ${params.ne}"
//log.info "         number of threads of execution                         : ${params.nthreads_SNP}"
//log.info "" }
log.info "         =================================================="
log.info "         =================================================="
log.info "         Output"
log.info "         =================================================="
log.info "         output FDR                      : ${params.output_FDR}"
if(params.SNPs && Gmodel){
log.info "         methQTL distance                : ${params.distance}" }
if(params.SNPs && GxE){
log.info "         number of kplots                : ${params.kplots}" }
log.info ""    
log.info "         =================================================="
log.info "         RUN NAME: ${workflow.runName}"
log.info ""



////////////////////
// STAGE CHANNELS //
////////////////////

// STAGE SAMPLES CHANNEL
samples_channel = Channel
    .from(file("${params.samples}").readLines())
    .ifEmpty{ exit 1, "ERROR: samples file is missing. Please remember to use the --samples parameter." }
    .map { line ->
        def field = line.toString().tokenize('\t').take(1)
        return tuple(field[0].replaceAll("\\s",""))}

//stage gff file with goa channel
//goa = "${params.goa}"
//species= "${params.species}"


// STAGE TEST PROFILE 
if ( workflow.profile.tokenize(",").contains("test") ){

    include {check_test_data} from './lib/functions.nf' params(CpGPaths: params.CpGPaths, CHGPaths: params.CHGPaths, CpGPaths_DMRs: params.CpGPaths_DMRs, CHGPaths_DMRs: params.CHGPaths_DMRs, SNPPaths: params.SNPPaths)
    (CpG, CHG, CpG_DMRs, CHG_DMRs, SNPs) = check_test_data(params.CpGPaths, params.CHGPaths, params.CpGPaths_DMRs, params.CHGPaths_DMRs, params.SNPPaths)

    CpG_DMPs = CpG_DMRs
    CHG_DMPs = CHG_DMRs
    
    CHH = Channel.empty()
    CHH_DMPs = Channel.empty()
    CHH_DMRs = Channel.empty()

} else {

    // STAGE INPUT CHANNELS
    CpG = params.noCpG  ? Channel.empty() :
        Channel
            .fromFilePairs( CpG_path, size: 1)
            .ifEmpty{ exit 1, "ERROR: No input found for CpG *.bedGraph files: ${params.input}/CpG\n\n \
            -Please check files exist or specify --noCpG \n \
            -Please check sample names match: ${samples}"}
    CHG = params.noCHG  ? Channel.empty() :
        Channel
            .fromFilePairs( CHG_path, size: 1)
            .ifEmpty{ exit 1, "ERROR: No input found for CHG *.bedGraph files: ${params.input}/CHG\n\n \
            -Please check files exist or specify --noCHG \n \
            -Please check sample names match: ${samples}"}
    CHH = params.noCHH  ? Channel.empty() :
        Channel
            .fromFilePairs( CHH_path, size: 1)
            .ifEmpty{ exit 1, "ERROR: No input found for CHH *.bedGraph files: ${params.input}/CHH\n\n \
            -Please check files exist or specify --noCHH \n \
            -Please check sample names match: ${samples}"}

    // STAGE DMP CHANNELS
    CpG_DMPs = params.noCpG  ? Channel.empty() : !params.DMPs ? Channel.empty() :
        Channel
            .fromFilePairs( CpG_path_DMPs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No DMP input found for CpG *.bed files: ${params.DMPs}/CpG\n\n \
            -Please check files exist or specify --noCpG"}
    CHG_DMPs = params.noCHG  ? Channel.empty() : !params.DMPs ? Channel.empty() :
        Channel
            .fromFilePairs( CHG_path_DMPs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No DMP input found for CHG *.bed files: ${params.DMPs}/CHG\n\n \
            -Please check files exist or specify --noCHG"}
    CHH_DMPs = params.noCHH  ? Channel.empty() : !params.DMPs ? Channel.empty() :
        Channel
            .fromFilePairs( CHH_path_DMPs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No DMP input found for CHH *.bed files: ${params.DMPs}/CHH\n\n \
            -Please check files exist or specify --noCHH"}

    // STAGE DMR CHANNELS
    CpG_DMRs = params.noCpG  ? Channel.empty() : !params.DMRs ? Channel.empty() :
        Channel
            .fromFilePairs( CpG_path_DMRs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No DMR input found for CpG *.bed files: ${params.DMRs}\n\n \
            -Please check files exist or specify --noCpG"}
    CHG_DMRs = params.noCHG  ? Channel.empty() : !params.DMRs ? Channel.empty() :
        Channel
            .fromFilePairs( CHG_path_DMRs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No DMR input found for CHG *.bed files: ${params.DMRs}\n\n \
            -Please check files exist or specify --noCHG"}
    CHH_DMRs = params.noCHH  ? Channel.empty() : !params.DMRs ? Channel.empty() :
        Channel
            .fromFilePairs( CHH_path_DMRs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No DMR input found for CHH *.bed files: ${params.DMRs}\n\n \
            -Please check files exist or specify --noCHH"}

    // STAGE SNPs CHANNEL
    globs = ["${params.SNPs}","${params.SNPs}/**.${params.extension}"]
    SNPs = !params.SNPs ? Channel.empty() : 
        Channel
            .fromFilePairs( globs, size: 1, type: "file" )
            .ifEmpty{ exit 1, "ERROR: No input found for SNP file(s): ${params.SNPs}\n\n \
            For single-sample vcfs:\n \
            -Please check files exist: ${params.SNPs}/vcf/*.${params.extension}\n \
            -Please check sample names match: ${samples}\n \
            -Please check given file extension: ${params.extension}"}
    
    
    //GOA = file("${params.GOA}", checkIfExists: true, glob: false)
    //species = file("${params.species}", checkIfExists: true, glob: false)

      
      
      
}

// METHYLATION CALLS
CpG_single = CpG.combine(samples_channel, by: 0).map{tuple("CpG", "bedGraph", *it)}
CHG_single = CHG.combine(samples_channel, by: 0).map{tuple("CHG", "bedGraph", *it)}
CHH_single = CHH.combine(samples_channel, by: 0).map{tuple("CHH", "bedGraph", *it)}
single_channel = CpG_single.mix(CHG_single,CHH_single)

// METHYLATION DMPs
CpG_DMPs_single = CpG_DMPs.map{tuple("CpG", "DMPs", *it)}
CHG_DMPs_single = CHG_DMPs.map{tuple("CHG", "DMPs", *it)}
CHH_DMPs_single = CHH_DMPs.map{tuple("CHH", "DMPs", *it)}
DMPs_channel = CpG_DMPs_single.mix(CHG_DMPs_single,CHH_DMPs_single)

// METHYLATION DMRs
CpG_DMRs_single = CpG_DMRs.map{tuple("CpG", "DMRs", *it)}
CHG_DMRs_single = CHG_DMRs.map{tuple("CHG", "DMRs", *it)}
CHH_DMRs_single = CHH_DMRs.map{tuple("CHH", "DMRs", *it)}
DMRs_channel = CpG_DMRs_single.mix(CHG_DMRs_single,CHH_DMRs_single)

// STAGE FINAL INPUTS
samples = file("${params.samples}", checkIfExists: true)
input_channel = single_channel.mix(DMPs_channel, DMRs_channel)

////////////////////
// BEGIN PIPELINE //
////////////////////

// INCLUDES
//include {parsing;split_scaffolds;calculate_FDR;qqPlot;GO_analysis} from './lib/process.nf' params(params)
include {parsing;split_scaffolds;calculate_FDR;qqPlot} from './lib/process.nf' params(params)
include {filtering;bedtools_unionbedg;bedtools_filtering;bedtools_sorting;bedtools_intersect;filter_regions;bedtools_merge;average_over_regions;GEM_Emodel;manhattan;} from './lib/GEM_Emodel.nf' params(params)
include {tabix;bcftools;vcftools_missing;vcftools_extract;GEM_Gmodel;GEM_GxEmodel;dotPlot;topKplots} from './lib/GEM_Gmodel.nf' params(params)
include {checkLines} from './lib/functions.nf'

// SUB-WORKFLOWS
workflow 'EWAS' {

    // get the initial files / Channels
    take:
        samples
        input_channel
        SNPs
        //GOA
        //species
        //goa
      

    // outline workflow
    main:
        // parse the samples.tsv file to get cov.txt and env.txt
        parsing(samples)

        // perform filtering on individual files
        filtering(input_channel)
        // stage channels for downstream processes
        bedGraph_combined = filtering.out.filter{it[1] == "bedGraph"}.groupTuple()
        DMPs_combined = filtering.out.filter{it[1] == "DMPs"}.groupTuple()
        DMRs_combined = filtering.out.filter{it[1] == "DMRs"}.groupTuple()
        bedtools_input = bedGraph_combined.mix(DMPs_combined, DMRs_combined)

        // bedtools_unionbedg for taking the union set in each context
        bedtools_unionbedg(bedtools_input.filter{ it[3].size() > 1 })
        unfiltered = bedtools_unionbedg.out.filter{ it[1] == "bedGraph" }.map{ tuple(it[0], "unfilter", *it[2..-1]) }

        // filtering bedGraph union bed files based on SD and NA
        bedtools_filtering(bedGraph_combined.filter{ it[3].size() == 1 }.mix(bedtools_unionbedg.out.filter{ it[1] == "bedGraph" }))
        bedtools_filtering_output = bedtools_filtering.out.filter{ checkLines(it[3]) > 1 }
        bedtools_filtering.out.filter{ checkLines(it[3]) <= 1 }.subscribe {
            log.warn "bedtools_filtering: no data left to analyse after filtering: ${it[0]}.${it[1]}.bed"
        }

        // sorting on union bed files
        sort_DMPs = DMPs_combined.filter{ it[3].size() == 1 }.mix(bedtools_unionbedg.out.filter{ it[1] == "DMPs" })
        sort_DMRs = DMRs_combined.filter{ it[3].size() == 1 }.mix(bedtools_unionbedg.out.filter{ it[1] == "DMRs" })
        bedtools_sorting(bedtools_filtering_output.mix(sort_DMPs, sort_DMRs, unfiltered))
        
        // stage channels for downstream processes
        intersect_DMPs = bedtools_sorting.out.filter{ it[1] == "DMPs" }
        intersect_DMRs = bedtools_sorting.out.filter{ it[1] == "DMRs" }
        bedGraph_DMPs = bedtools_sorting.out.filter{it[1] == "bedGraph"}.combine(intersect_DMPs, by: 0)
        bedGraph_DMRs = bedtools_sorting.out.filter{it[1] == "bedGraph"}.combine(intersect_DMRs, by: 0)
        unfilter_DMRs = bedtools_sorting.out.filter{it[1] == "unfilter"}.combine(intersect_DMRs, by: 0)

        // bedtools_intersect for intersecting individual methylation info based on DMPs/DMRs
        bedtools_intersect(bedGraph_DMPs.mix(bedGraph_DMRs))

        // filter regions based on bootstrap values
        filter_regions(bedGraph_DMRs)
        filter_regions_output = filter_regions.out.filter{ checkLines(it[4]) > 1 }
        filter_regions.out.filter{ checkLines(it[4]) <= 1 }.subscribe {
            log.warn "filter_regions: no data left to analyse after filtering: ${it[0]}.${it[1]}.bed"
        }
        // bedtools_merge for optionally combining filtered sub-regions
        bedtools_merge(filter_regions_output)
        // average_over_regions for calculating average methylation over defined regions
        average_over_regions(filter_regions_output.mix(bedtools_merge.out))
        average_over_regions_output = average_over_regions.out.filter{ checkLines(it[2]) > 1 }
        average_over_regions.out.filter{ checkLines(it[2]) <= 1 }.subscribe {
            log.warn "average_over_regions: no data left to analyse after intersection: ${it[0]}.${it[1]}.bed"
        }

        // stage channels for downstream processes
        bedGraph_channel = params.all || (!params.DMPs && !params.DMRs) ? bedtools_sorting.out.filter{it[1] == "bedGraph"} : Channel.empty()
        meth_channel = bedGraph_channel.mix(bedtools_intersect.out, average_over_regions.out)
        split_scaffolds(meth_channel)
        meth_channel = bedGraph_channel.mix(bedtools_intersect.out, average_over_regions_output)

        // SNPs
        // index individual vcf files, optionally rename header
        tabix(SNPs)
        // merge files, normalise, validate sample names in header
        bcftools(samples, tabix.out[0].collect(), tabix.out[1].collect())
        // extract missing information
        vcftools_missing(bcftools.out)
        //SNP imputation with BEAGLE
        BEAGLE_SNP_Imputation(bcftools.out)
        // extract snps.txt for GEM_GModel
        vcftools_extract(bcftools.out)

        // run GEM on selected combination of inputs
        GEM_Emodel(split_scaffolds.out.transpose(), parsing.out[0], parsing.out[1])
        GEM_Gmodel(split_scaffolds.out.transpose(), vcftools_extract.out, parsing.out[1])
        GEM_GxEmodel(split_scaffolds.out.transpose(), vcftools_extract.out, parsing.out[2])
        
        // calculate FDR
        Emodel_txt = GEM_Emodel.out[0].collectFile().map{ tuple(it.baseName, it) }
        Emodel_log = GEM_Emodel.out[1].collectFile().map{ tuple(it.baseName, it) }
        Emodel_channel = Emodel_txt.combine(Emodel_log, by: 0)
            .map { it -> 
                List key = it[0].tokenize(".")
                context = key.init().join(".")
                type = key.last()
                return tuple("Emodel", it[0], context, type, it[1], it[2])
            }

        Gmodel_txt = GEM_Gmodel.out[0].collectFile().map{ tuple(it.baseName, it) }
        Gmodel_log = GEM_Gmodel.out[1].collectFile().map{ tuple(it.baseName, it) }
        Gmodel_channel = Gmodel_txt.combine(Gmodel_log, by: 0)
            .map { it -> 
                List key = it[0].tokenize(".")
                context = key.init().join(".")
                type = key.last()
                return tuple("Gmodel", it[0], context, type, it[1], it[2])
            }

        GxE_txt = GEM_GxEmodel.out[0].collectFile().map{ tuple(it.baseName, it) }
        GxE_log = GEM_GxEmodel.out[1].collectFile().map{ tuple(it.baseName, it) }
        GxE_channel = GxE_txt.combine(GxE_log, by: 0)
            .map { it -> 
                List key = it[0].tokenize(".")
                context = key.init().join(".")
                type = key.last()
                return tuple("GxE", it[0], context, type, it[1], it[2])
            }

        // eg. [Emodel, context, type, context.txt, [scaffold.txt, ...], [scaffold.log], ...]]
        calculate_FDR(Emodel_channel.mix(Gmodel_channel, GxE_channel))
        //calculate_FDR_output = calculate_FDR.out.filter{ checkLines(it[3]) > 1 }
        //calculate_FDR.out.filter{ checkLines(it[3]) <= 1 }.subscribe {
        //    log.warn "calculate_FDR: no data left to analyse after filtering: ${it[1]}.${it[2]}.txt"
        //}
        
        //GOA
        //GO_analysis(calculate_FDR.out)
        
        // visualisation
        qqPlot(calculate_FDR.out)
        //GO_analysis(goa, species, calculate_FDR.out[0].collect())
        //GO_analysis(goa)
        manhattan(calculate_FDR.out.filter{ it[0] == "Emodel" })
        dotPlot(calculate_FDR.out.filter{ it[0] == "Gmodel" })
        //kplots_channel = calculate_FDR.out.filter{ it[0] == "GxE" }.map{ it.tail() }.combine(GEM_GxEmodel.out[1].map{ tuple( it[0] + "." + it[1], it.last()) }.groupTuple(), by:0)

        GxE_plot = GEM_GxEmodel.out[2].collectFile().map{ tuple(it.baseName, it) }
        GxE_head = GEM_GxEmodel.out[3].map{ tuple( it[0] + "." + it[1], it[2]) }.unique{ it[0] }
        kplots_channel = calculate_FDR.out.filter{ it[0] == "GxE" }.map{ it.tail() }.combine(GxE_plot, by:0).combine(GxE_head, by:0)
        topKplots(kplots_channel, vcftools_extract.out, parsing.out[2])

/*
    // emit results
    emit:
        parsing_env = parsing.out[0]
        parsing_cov = parsing.out[1]
        parsing_gxe = GxE ? parsing.out[2] : Channel.empty()
        vcftools_extract_out = vcftools_extract.out
        BEAGLE_SNP_Imputation_out = BEAGLE_SNP_Imputation.out 

        bedtools_unionbedg_out = bedtools_sorting.out
        bedtools_intersect_out = bedtools_intersect.out
        average_over_regions_out = average_over_regions.out

        calculate_FDR_reg = calculate_FDR.out[0].filter{ it[2] == "region" || it[2] == "merged" }
        calculate_FDR_pos = calculate_FDR.out[0].filter{ it[2] != "region" && it[2] != "merged" }
        
        //GO_analysis_reg  = GO_analysis.out[0].filter{ it[2] == "region" || it[2] == "merged" }
        //GO_analysis_pos  = GO_analysis.out[0].filter{ it[2] == "region" || it[2] == "merged" }
        
        qqPlot_png_reg = qqPlot.out.filter{ it[0] == "region" || it[0] == "merged" }
        qqPlot_png_pos = qqPlot.out.filter{ it[0] != "region" && it[0] != "merged" }
        manhattan_png_reg = manhattan.out[0].filter{ it[0] == "region" || it[0] == "merged" }
        manhattan_png_pos = manhattan.out[0].filter{ it[0] != "region" && it[0] != "merged" }
        manhattan_zip_reg = manhattan.out[1].filter{ it[0] == "region" || it[0] == "merged" }
        manhattan_zip_pos = manhattan.out[1].filter{ it[0] != "region" && it[0] != "merged" }
        dotPlot_png_reg = dotPlot.out[0].filter{ it[0] == "region" || it[0] == "merged" }
        dotPlot_png_pos = dotPlot.out[0].filter{ it[0] != "region" && it[0] != "merged" }
        dotPlot_zip_reg = dotPlot.out[1].filter{ it[0] == "region" || it[0] == "merged" }
        dotPlot_zip_pos = dotPlot.out[1].filter{ it[0] != "region" && it[0] != "merged" }
        topKplots_png_reg = topKplots.out.filter{ it[0] == "region" || it[0] == "merged" }
        topKplots_png_pos = topKplots.out.filter{ it[0] != "region" && it[0] != "merged" }
*/

}

// MAIN WORKFLOW 
workflow {

    // call sub-workflows
    main:
        EWAS(samples, input_channel, SNPs)
/*
    // publish files
    publish:
        EWAS.out.parsing_env to: "${params.output}/input", mode: 'copy'
        EWAS.out.parsing_cov to: "${params.output}/input", mode: 'copy'
        EWAS.out.parsing_gxe to: "${params.output}/input", mode: 'copy'
        EWAS.out.vcftools_extract_out to: "${params.output}/input", mode: 'copy'

        EWAS.out.bedtools_unionbedg_out to: "${params.output}/input", mode: 'copy'
        EWAS.out.bedtools_intersect_out to: "${params.output}/positions", mode: 'copy'
        EWAS.out.average_over_regions_out to: "${params.output}/regions", mode: 'copy'

        EWAS.out.calculate_FDR_reg to: "${params.output}/regions", mode: 'copy'
        EWAS.out.calculate_FDR_pos to: "${params.output}/positions", mode: 'copy'
        
        //EWAS.out.GO_analysis_reg to: "${params.output}/regions/GOA", mode: 'copy'
        //EWAS.out.GO_analysis_pos to: "${params.output}/positions/GOA", mode: 'copy'

        EWAS.out.qqPlot_png_reg to: "${params.output}/regions", mode: 'copy'
        EWAS.out.qqPlot_png_pos to: "${params.output}/positions", mode: 'copy'
        EWAS.out.manhattan_png_reg to: "${params.output}/regions/Emodel", mode: 'copy'
        EWAS.out.manhattan_png_pos to: "${params.output}/positions/Emodel", mode: 'copy'
        EWAS.out.manhattan_zip_reg to: "${params.output}/regions/Emodel", mode: 'copy'
        EWAS.out.manhattan_zip_pos to: "${params.output}/positions/Emodel", mode: 'copy'
        EWAS.out.dotPlot_png_reg to: "${params.output}/regions", mode: 'copy'
        EWAS.out.dotPlot_png_pos to: "${params.output}/positions", mode: 'copy'
        EWAS.out.dotPlot_zip_reg to: "${params.output}/regions", mode: 'copy'
        EWAS.out.dotPlot_zip_pos to: "${params.output}/positions", mode: 'copy'
        EWAS.out.topKplots_png_reg to: "${params.output}/regions", mode: 'copy'
        EWAS.out.topKplots_png_pos to: "${params.output}/positions", mode: 'copy'
*/
}



//////////////////
// END PIPELINE //
//////////////////

// WORKFLOW TRACING # what to display when the pipeline finishes
// eg. with errors
workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

// eg. in general
workflow.onComplete {

    log.info ""
    log.info "         Pipeline execution summary"
    log.info "         ---------------------------"
    log.info "         Name         : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}"
    log.info "         Profile      : ${workflow.profile}"
    log.info "         Launch dir   : ${workflow.launchDir}"    
    log.info "         Work dir     : ${workflow.workDir} ${workflow.success && !params.debug ? "(cleared)" : ""}"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""

    // run a small clean-up script to remove "work" directory after successful completion 
    if (workflow.success && !params.debug) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}
