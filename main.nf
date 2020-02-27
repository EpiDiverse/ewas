#!/usr/bin/env nextflow

// DSL2 BRANCH
nextflow.preview.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\

         =================================================
          E P I D I V E R S E - E W A S   P I P E L I N E
         =================================================
         ~ version ${workflow.manifest.version}

         Usage: 
              nextflow run epidiverse/template [OPTIONS]...

         Options: GENERAL
              --input [path/to/input/dir]     [REQUIRED] Write a line about the input file(s) format

              --reference [path/to/ref.fa]    [REQUIRED] Write a line aboue the reference file(s) format

              --output [STR]                  A string that can be given to name the output directory [default: output]


         Options: MODIFIERS
              --process1                      A boolean (true/false) parameter which alters the default behaviour [default: off]

              --process2                      Another boolean parameter. Note: default boolean params must always be false [default: off]

              --paired                        Enable paired-end mode. Another boolean parameter. [default: off]

         Options: ADDITIONAL
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    


         Example: 
              nextflow run epidiverse/template \
              --input /path/to/input/dir/or/file \
              --reference /path/to/ref.fa
              --output output

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


// DEFINE PATHS
CpG_path = "${params.input}/*/bedGraph/*_CpG.bedGraph"
CHG_path = "${params.input}/*/bedGraph/*_CHG.bedGraph"
CHH_path = "${params.input}/*/bedGraph/*_CHH.bedGraph"
  
CpG_path_DMPs = "${params.DMPs}/CpG/metilene/*" 
CHG_path_DMPs = "${params.DMPs}/CHG/metilene/*"
CHH_path_DMPs = "${params.DMPs}/CHH/metilene/*"

CpG_path_DMRs = "${params.DMRs}/CpG/metilene/*" 
CHG_path_DMRs = "${params.DMRs}/CHG/metilene/*"
CHH_path_DMRs = "${params.DMRs}/CHH/metilene/*"

//snp_path_mult = "${params.SNPs}"
//snp_path = "${params.SNPs}/vcf/*.vcf"


// PARAMETER CHECKS
if( !params.input && (params.DMPs || params.DMRs) ){error "ERROR: "}
if( params.noCpG && params.noCHG && params.noCHH ){error "ERROR: please specify at least one methylation context for analysis"}


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
log.info "         input dir    : ${params.input}"
log.info "         output dir   : ${params.output}"
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
    .ifEmpty{ exit 1, "ERROR: samples file is missing or invalid. Please remember to use the --samples parameter." }
    .map { line ->
        def field = line.toString().tokenize('\t').take(1)
        return tuple(field[0].replaceAll("\\s",""))}


// STAGE BAM FILES FROM TEST PROFILE # this establishes the test data to use with -profile test
if ( workflow.profile.tokenize(",").contains("test") ){

        include check_test_data from './libs/functions.nf' params(BAMPaths: params.BAMPaths)
        BAM = check_test_data(params.BAMPaths)

} else {

    // STAGE INPUT CHANNELS
    CpG = params.noCpG  ? Channel.empty() : !params.input ? Channel.empty() :
        Channel
            .fromFilePairs( CpG_path, size: 1)
            .ifEmpty{ exit 1, "ERROR: No input found for CpG files: ${params.input}\n"}
    CHG = params.noCHG  ? Channel.empty() : !params.input ? Channel.empty() :
        Channel
            .fromFilePairs( CHG_path, size: 1)
            .ifEmpty{ exit 1, "ERROR: No input found for CHG files: ${params.input}\n"}
    CHH = params.noCHH  ? Channel.empty() : !params.input ? Channel.empty() :
        Channel
            .fromFilePairs( CHH_path, size: 1)
            .ifEmpty{ exit 1, "ERROR: No input found for CHH files: ${params.input}\n"}

    // STAGE DMP CHANNELS
    CpG_DMPs = params.noCpG  ? Channel.empty() : !params.DMPs ? Channel.empty() :
        Channel
            .fromFilePairs( CpG_path_DMPs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No input found for CpG files: ${params.DMPs}\n"}
    CHG_DMPs = params.noCHG  ? Channel.empty() : !params.DMPs ? Channel.empty() :
        Channel
            .fromFilePairs( CHG_path_DMPs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No input found for CHG files: ${params.DMPs}\n"}
    CHH_DMPs = params.noCHH  ? Channel.empty() : !params.DMPs ? Channel.empty() :
        Channel
            .fromFilePairs( CHH_path_DMPs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No input found for CHH files: ${params.DMPs}\n"}

    // STAGE DMR CHANNELS
    CpG_DMRs = params.noCpG  ? Channel.empty() : !params.DMRs ? Channel.empty() :
        Channel
            .fromFilePairs( CpG_path_DMRs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No input found for CpG files: ${params.DMRs}\n"}
    CHG_DMRs = params.noCHG  ? Channel.empty() : !params.DMRs ? Channel.empty() :
        Channel
            .fromFilePairs( CHG_path_DMRs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No input found for CHG files: ${params.DMRs}\n"}
    CHH_DMRs = params.noCHH  ? Channel.empty() : !params.DMRs ? Channel.empty() :
        Channel
            .fromFilePairs( CHH_path_DMRs, size: 1, type: "dir")
            .ifEmpty{ exit 1, "ERROR: No input found for CHH files: ${params.DMRs}\n"}

}

// METHYLATION CALLS
CpG_single = CpG.combine(samples_channel).map{tuple("CpG", "bedGraph", *it)}
CHG_single = CHG.combine(samples_channel).map{tuple("CHG", "bedGraph", *it)}
CHH_single = CHH.combine(samples_channel).map{tuple("CHH", "bedGraph", *it)}
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

// INCLUDES # here you must give the relevant process files from the libs directory 
include './libs/process.nf' params(params)

// SUB-WORKFLOWS
workflow 'EWAS' {

    // get the initial files / Channels
    get:
        samples
        input_channel

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

        // bedtools_unionbedg for taking the union set in each context
        bedtools_unionbedg(bedGraph_combined.mix(DMPs_combined, DMRs_combined))

        // stage channels for downstream processes
        bedGraph_DMPs = bedtools_unionbedg.out.filter{it[1] == "bedGraph"}.combine(bedtools_unionbedg.out.filter{it[1] == "DMPs"}, by: 0)
        bedGraph_DMRs = bedtools_unionbedg.out.filter{it[1] == "bedGraph"}.combine(bedtools_unionbedg.out.filter{it[1] == "DMRs"}, by: 0)
        intersect_channel = bedGraph_DMPs.mix(bedGraph_DMRs)

        // bedtools_intersect for intersecting individual methylation info based on DMPs/DMRs
        bedtools_intersect(bedGraph_DMPs.mix(bedGraph_DMRs))

        // bedtools_merge for optionally combining filtered sub-regions
        bedtools_merge(bedGraph_DMRs)

        // average_over_regions for calculating average methylation over defined regions
        average_over_regions(bedGraph_DMRs.mix(bedtools_merge.out))

        // run GEM_Emodel on selected combination of inputs
        emodel_channel = bedtools_unionbedg.out.filter{it[1] == "bedGraph"}.mix(bedtools_intersect.out, average_over_regions.out)
        GEM_Emodel(emodel_channel, parsing.out[0], parsing.out[1])

        // additional steps for further downstream processing 
        // ..

    // Emit results
    emit:
        parsing_env = parsing.out[0]
        parsing_cov = parsing.out[1]
        bedtools_unionbedg_out = bedtools_unionbedg.out
        bedtools_intersect_out = bedtools_intersect.out
        average_over_regions_out = average_over_regions.out
        //GEM_Emodel_filtered_pos = GEM_Emodel.out[0].filter{it[1] != "region"}
        //GEM_Emodel_filtered_reg = GEM_Emodel.out[0].filter{it[1] == "region"}
        //GEM_Emodel_full_pos = GEM_Emodel.out[1].filter{it[1] != "region"}
        //GEM_Emodel_full_reg = GEM_Emodel.out[1].filter{it[1] == "region"}
        //GEM_Emodel_jpg_pos = GEM_Emodel.out[2].filter{it[1] != "region"}
        //GEM_Emodel_jpg_reg = GEM_Emodel.out[2].filter{it[1] == "region"}
}

// MAIN WORKFLOW 
workflow {

    // call sub-workflows eg. WORKFLOW(Channel1, Channel2, Channel3, etc.)
    main:
        EWAS(samples, input_channel)

    // define how the emitted files from each sub-workflow should be "published" for the user to see
    publish:
        EWAS.out.parsing_env to: "${params.output}", mode: 'copy'
        EWAS.out.parsing_cov to: "${params.output}", mode: 'copy'
        EWAS.out.bedtools_unionbedg_out to: "${params.output}/input", mode: 'copy'
        EWAS.out.bedtools_intersect_out to: "${params.output}/positions", mode: 'copy'
        EWAS.out.average_over_bed_out to: "${params.output}/regions", mode: 'copy'
        //EWAS.out.GEM_Emodel_filtered_pos to: "${params.output}/positions/Emodel", mode: 'copy'
        //EWAS.out.GEM_Emodel_filtered_reg to: "${params.output}/regions/Emodel", mode: 'copy'
        //EWAS.out.GEM_Emodel_full_pos to: "${params.output}/positions/Emodel", mode: 'copy'
        //EWAS.out.GEM_Emodel_full_reg to: "${params.output}/regions/Emodel", mode: 'copy'
        //EWAS.out.GEM_Emodel_jpg_pos to: "${params.output}/positions/Emodel", mode: 'copy'
        //EWAS.out.GEM_Emodel_jpg_reg to: "${params.output}/regions/Emodel", mode: 'copy'

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
    log.info "         Work dir     : ${workflow.workDir} ${params.debug ? "" : "(cleared)" }"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""

    // run a small clean-up script to remove "work" directory after successful completion 
    if (params.debug == false && workflow.success) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}