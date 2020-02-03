#!/usr/bin/env nextflow

// DSL2 BRANCH
nextflow.preview.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\

         ===============================================
          E P I D I V E R S E - X X X   P I P E L I N E
         ===============================================
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
         ===============================================
          E P I D I V E R S E - X X X   P I P E L I N E
         ===============================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}


// DEFINE PATHS # these are strings which are used to define input Channels,
// but they are specified here as they may be referenced in LOGGING
bam_path = "${params.input}/*/*.bam"

// SET CONDITIONALS # some user parameters may influence each other
// eg. here the default behaviour is "process1" but this is switched off
// when the user specifies --process2. Alternatively, the user may specify
// both --process1 and --process2 to perform both.
if( !params.process2 ){
    process1 = true
    process2 = false
} else if( params.process1 ){
    process1 = true
    process2 = true
} else {
    process1 = false
    process2 = true
}

// MORE CONDITIONALS # eg. here the input reference is only checked in 
// "process1" mode
if( process1 ){
    fasta = file("${params.reference}", checkIfExists: true, glob: false)
    fai = file("${params.reference}.fai", checkIfExists: true, glob: false)
} else {
    fasta = false
    fai = false
}

// PRINT STANDARD LOGGING INFO
log.info ""
log.info "         ================================================"
log.info "          E P I D I V E R S E - S N P    P I P E L I N E"
if(params.debug){
log.info "         (debug mode enabled)"
log.info "         ================================================" }
else {
log.info "         ================================================" }
log.info "         ~ version ${workflow.manifest.version}"
log.info ""
log.info "         input dir    : ${params.input}"
log.info "         reference    : ${params.reference ? "${params.reference}" : "-"}"
log.info "         output dir   : ${params.output}"
log.info "         process1     : ${process1 ? "enabled" : "disabled"}"
log.info "         process2     : ${process2 ? "enabled" : "disabled"}"
log.info ""
log.info "         ================================================"
log.info "         RUN NAME: ${workflow.runName}"
log.info ""



////////////////////
// STAGE CHANNELS //
////////////////////

/*
 *   Channels are where you define the input for the different
 *    processes which make up a pipeline. Channels indicate
 *    the flow of data, aka the "route" that a file will take.
 */

// STAGE BAM FILES FROM TEST PROFILE # this establishes the test data to use with -profile test
if ( workflow.profile.tokenize(",").contains("test") ){

        include check_test_data from './libs/functions.nf' params(BAMPaths: params.BAMPaths)
        BAM = check_test_data(params.BAMPaths)

} else {

    // STAGE BAM CHANNELS # this defines the normal input when test profile is not in use
    BAM = Channel.fromPath(bam_path)
        .ifEmpty{ exit 1, "ERROR: cannot find valid *.bam files in dir: ${params.input}\n"}
        .map{ tuple(it.parent.name, it) }
        .take(params.take.toInteger())
}

// ERROR HANDLING # for Channels which may contain bad data
// eg. here the pipeline stops in "process2" mode with fewer than three samples
if( process2 ){

    BAM
        .count()
        .subscribe{int c ->
            if( c <= 2 ){
                error "ERROR: process2 is only possible with a minimum of three samples"
                exit 1
            }
        }

}


////////////////////
// BEGIN PIPELINE //
////////////////////

/*
 *   Workflows are where you define how different processes link together. They
 *    may be modularised into "sub-workflows" which must be named eg. 'TEMPLATE'
 *    and there must always be one MAIN workflow to link them together, which
 *    is always unnamed.
 */


// INCLUDES # here you must give the relevant process files from the libs directory 
include './libs/process.nf' params(params)

// SUB-WORKFLOWS
workflow 'TEMPLATE' {

    // get the initial Channels
    get:
        BAM
        fasta
        fai

    // define how the different processes lead into each other with
    // eg. process(input1, input2, etc.)
    // eg. process.out[0], process.out[1], etc.
    main:
        // eg. samtools sort + index
        preprocessing(BAM)

        // process outputs can be "forked" into multiple Channels for different processes 
        process1_channel = process1 ? preprocessing.out.map{tuple("process1", *it)} : Channel.empty()
        process2_channel = process2 ? preprocessing.out.map{tuple("process2", *it)} : Channel.empty()

        // now process1 will only run with process1_channel, process2 with process2_channel
        process1(process1_channel, fasta, fai)
        process2(process2_channel)

        // ... they can be combined back together again in process3. You can do whatever you like.
        process3(process1.out.mix(process2.out))


    // Emit the relevant output which should be "published" for the user to see
    // index numbers 0,1,etc. refer to different outputs defined for processes in process.nf
    emit:
        process1_output = process1.out
        process2_output_bam = process2.out[0]
        process2_output_txt = process2.out[1]
        process3_output = process3.out
}

// MAIN WORKFLOW 
workflow {

    // call sub-workflows eg. WORKFLOW(Channel1, Channel2, Channel3, etc.)
    main:
        TEMPLATE(BAM, fasta, fai)

    // define how the emitted files from each sub-workflow should be "published" for the user to see
    publish:
        TEMPLATE.out.process1_output to: "${params.output}/bam", mode: 'copy'
        TEMPLATE.out.process2_out_bam to: "${params.output}/bam", mode: 'copy'
        TEMPLATE.out.process2_out_txt to: "${params.output}", mode: 'move'
        TEMPLATE.out.process3_out to: "${params.output}/fastq", mode: 'copy'
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