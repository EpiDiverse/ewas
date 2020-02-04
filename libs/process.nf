#!/usr/bin/env nextflow

if (params.DMP && params.DMR) {
    println 'ERROR: Please specify at most one input type!'
}
else if (params.DMP) {   
    CpG_metilene_path        = "${params.DMP}/CpG/metilene/*" 
    CHG_metilene_path        = "${params.DMP}/CHG/metilene/*"
    CHH_metilene_path        = "${params.DMP}/CHH/metilene/*"   
}
else if (params.DMR) {   
    CpG_metilene_path        = "${params.DMR}/CpG/metilene/*" 
    CHG_metilene_path        = "${params.DMR}/CHG/metilene/*"
    CHH_metilene_path        = "${params.DMR}/CHH/metilene/*"   
}


// DEFINE MODEL TYPE FOR GEM RUN
def modelLine = ""
if( params.GEM_Emodel == true ){ modelLine += "GEM_Emodel," }
if( params.GEM_Gmodel == true ){ modelLine += "GEM_Gmodel," }
if( params.GEM_GXEmodel == true ){ modelLine += "GEM_GXEmodel," }
if( (params.GEM_Emodel == false) && (params.GEM_Gmodel == false) && (params.GEM_GXEmodel == false) ){error "ERROR: please specify at least one model type for GEM run"}



// DEFINE CONTEXTLINE FOR INPUT PATH
def contextLine = ""
if( params.noCpG == false ){ contextLine += "CpG," }
if( params.noCHG == false ){ contextLine += "CHG," }
if( params.noCHH == false ){ contextLine += "CHH," }
if( params.noCpG && params.noCHG && params.noCHH ){error "ERROR: please specify at least one methylation context for analysis"}