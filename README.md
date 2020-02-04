EpiDiverse-EWAS Pipeline
========================

This pipeline is configured for running on the EPI server cluster.

You must first make the pipeline available using this command:

	nextflow pull https://bitbucket.org/epidiverse/ewas

GEM is an R tool suite for performing epigenome wide association studies (EWAS). GEM provides three major functions named GEM_Emodel, 

GEM_Gmodel and GEM_GxEmodel to study the interaction of Gene, Environment and Methylation (GEM).



GEM_Emodel is used to find the association between methylation and environmental factor genome widely form using the following command:

Methylation calls:

	nextflow run epidiverse/ewas --input [path/to/bedGraph/files] --samples [path/to/samples.tsv] --GEM_Emodel --profile epi,custom

DMPs:

    nextflow run epidiverse/ewas --input [path/to/bedGraph/files] --samples [path/to/samples.tsv] --GEM_Emodel --DMP [path/to/DMPs]
    --profile epi,custom

DMRs:

    nextflow run epidiverse/ewas --input [path/to/bedGraph/files] --samples [path/to/samples.tsv] --GEM_Emodel --DMR [path/to/DMRs]
    --profile epi,custom



GEM_Gmodel is to create a methQTL genome-wide map using the following command:

	NOT IMPLEMENTED YET


GEM_GxEmodel is to test ability of the interaction of gene and environmental factor to predict DNA methylation level using the following command:

	NOT IMPLEMENTED YET


Example of a sample.tsv file:

    #sample_identifier             env   cov
    PI_XX_NN_NN_X_YYMMDD_X_1		32	  1
    PI_XX_NN_NN_X_YYMMDD_X_2		21	  1
    PI_XX_NN_NN_X_YYMMDD_X_3		10	  1







Pipeline Parameters
-------------------

	 INPUT / OUTPUT

  
              --input                         [REQUIRED] Specify the path to the directory containing each sample output from the WGBS 
                                              pipeline to be taken forward for analysis. All the subdirectories must correspond to sample 
                                              names in the provided samples.tsv file, and contain within them a bedGraph directories with 
                                              files in '*.bedGraph' format. 
                                              
                                              
              --DMP                           Specify the path to the directory containing each sample output from DMR pipeline to run EWAS 
                                              analysis with DMPs. [default: off]
                                             
         
              --DMR                           Specify the path to the directory containing each sample output from DMR pipeline to run EWAS 
                                              analysis with DMRs. [default: off]
              
              
              --samples [path/to/samples.tsv] [REQUIRED] Specify the path to the "samples.tsv" file containing information 
                                              regarding sample names and corresponding groupings/replicates. The file must contain
                                              three tab-separated columns: 1) sample names, corresponding to subdirectories in the
                                              --input directory, 2) group names, 3) replicates names, 4) environment values and 
                                              5) covariate values

              
              --snp                           [path/to/snp_file] [REQUIRED for GEM_G and GEM_GXE Models run]. A subset with genotype data 
                                              encoded as 1,2,3 for major allele homozygote (AA), heterozygote (AB) and minor allele homozygote 
                                              (BB) for all SNPs across all samples.  

              --output [STR]                  A string that will be used as the name for the output results directory, which will be generated 
                                              in the working directory. This directory will contain sub-directories for each set of reads analysed 
                                              during the pipeline. [default: 'ewas']

              --GEM_Gmodel                   Specify this parameter to run GEM with Gmodel. GEM_Gmodel creates a methQTL genome-wide map. 
                                             [default: off]
                                             
                                             
              --GEM_Emodel                   Specify this parameter to run GEM with Emodel. GEM_Emodel  aims to find the association between 
                                             methylation and environmental factor genome widely.  [default: off]
              
              --GEM_GXEmodel                 Specify this parameter to run GEM with GXEmodel. GEM_GXEmodel tests the ability of the interaction 
                                             of gene and environmental factor to predict DNA methylation level.  [default: off] 

                                             
      PIPELINE OPTIONS                                      
											 
              
              --noCpG                        Disables DMR analysis in CpG context. [default: off] 
              --noCHG                        Disables DMR analysis in CHG context. [default: off]
              --noCHH                        Disables DMR analysis in CHH context. [default: off]
              

    
	
	FILTERING
    
	
	          --sig                         Specify the maximum p-value threshold for filtering EWAS post-analysis. [default: 0.05]
              
              --FDR                         Specify the maximum FDR threshold for filtering EWAS post-analysis. [default: 0.05]

        
              --cov                         Specify the minimum coverage threshold to filter methylated positions before running
                                            the EWAS analyses with all input types. [default: 0] 
   
   
    DEVELOPER USE
        
	          --debug                       Prevent nextflow from clearing the cache on completion of the pipeline. This includes '.work' 
                                            and '.nextflow' directories and any log files that have been created. [default: off]
                                            
Pipeline Overview
-----------------

A general overview of the pipeline is given below. 


Output Directory Structure
--------------------------

An overview of the output directory structure produced by the pipeline. 

![picture](workflow/EWAS_Output_dir_second_draft.jpg)

Input Directory Structure
--------------------------

An overview of the input directory structure produced by the pipeline. 

![picture](workflow/EWAS_input_dir_second_draft.jpg)



Bitbucket Credentials (DEVELOPER USE)
-------------------------------------

In your `$HOME` directory, you can create a file to inform nextflow of your bitbucket credentials in order to
avoid typing `-user user@e-mail.com` and your password during the pipeline execution. A `template_scm.txt` file
has been included in the repository, where you can save your login credentials and store it under
`$HOME/.nextflow/scm` on the server cluster.