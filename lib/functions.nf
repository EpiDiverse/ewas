#!/usr/bin/env nextflow
// This file for loading custom functions into the main.nf script (separated for portability)

// FUNCTION TO LOAD DATASETS IN TEST PROFILE
def check_test_data(CpGPaths, CHGPaths, CpGPaths_DMPs, CHGPaths_DMPs, CpGPaths_DMRs, CHGPaths_DMRs, SNPPaths) {

    // STAGE INPUT CHANNELS
    CpG = Channel.from(CpGPaths)
        .map { row -> [ row[0], file(row[1]) ] }
        .ifEmpty { exit 1, "test profile CpGPaths was empty - no input files supplied" }

    CHG = Channel.from(CHGPaths)
        .map { row -> [ row[0], file(row[1]) ] }
        .ifEmpty { exit 1, "test profile CHGPaths was empty - no input files supplied" }

    // STAGE DMR CHANNELS
    CpG_DMRs = Channel.from(CpGPaths_DMRs)
        .map { row -> [ row[0], file(row[1]) ] }
        .ifEmpty { exit 1, "test profile CpGPaths_DMRs was empty - no input files supplied" }

    CHG_DMRs = Channel.from(CHGPaths_DMRs)
        .map { row -> [ row[0], file(row[1]) ] }
        .ifEmpty { exit 1, "test profile CHGPaths_DMRs was empty - no input files supplied" }

    // STAGE DMP CHANNELS
    CpG_DMPs = Channel.from(CpGPaths_DMPs)
        .map { row -> [ row[0], file(row[1]) ] }
        .ifEmpty { exit 1, "test profile CpGPaths_DMPs was empty - no input files supplied" }

    CHG_DMPs = Channel.from(CHGPaths_DMPs)
        .map { row -> [ row[0], file(row[1]) ] }
        .ifEmpty { exit 1, "test profile CHGPaths_DMPs was empty - no input files supplied" }
        
    // STAGE SNPs CHANNEL
    SNPs = Channel.from(SNPPaths)
        .map { row -> [ row[0], file(row[1]) ] }
        .ifEmpty { exit 1, "test profile SNPPaths was empty - no input files supplied" }

    // Return channels
    return tuple(CpG, CHG, CpG_DMPs, CHG_DMPs, CpG_DMRs, CHG_DMRs, SNPs)
}


def checkLines(myFile) {
    int count = 0
    String line = null
    myFile.withReader {
        while( count <= 1 && ( ( line = it.readLine() ) != null )) {count++}
    }
    return count
}
