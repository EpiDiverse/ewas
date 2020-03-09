# EpiDiverse-EWAS Usage
This document describes the parameter options used by the pipeline.

* [Running the pipeline](#running-the-pipeline)
* [Inputs and outputs](#inputs-and-outputs)
    * [`--input`](#--input-arg-required)
    * [`--samples`](#--samples-arg-required)
    * [`--DMPs`](#--dmps-arg)
    * [`--DMRs`](#--dmrs-arg)
    * [`--SNPs`](#--snps-arg)
    * [`--extension`](#--extension-arg)  
    * [`--output`](#--output-arg)
* [Model Decision](#model-decision)
    * [`--Emodel`](#--mmodel)
    * [`--Gmodel`](#--gmodel)
    * [`--GxE`](#--gxe) 
    * [`--noCpG`](#--nocpg)
    * [`--noCHG`](#--nochg)
    * [`--noCHH`](#--nochh)
* [Inputs Filtering](#input-filtering)
    * [`--coverage`](#--coverage-arg)
    * [`--input_FDR`](#--input-fdr-arg)
    * [`--proportion`](#--proportion-arg)
    * [`--merge`](#--merge)
* [SNP Filtering](#--SNP-Filtering) 
    * [`--max_missing`](#--max-missing)
    * [`--mac`](#--mac)
    * [`--minQ`](#--minQ)
* [Output Filtering](#--output-filtering)
    * [`--output_FDR`](#--output-fdr-arg)
    * [`--Emodel_pv`](#--emodel-pv-arg)
    * [`--Gmodel_pv`](#--gmodel-pv-arg)
    * [`--GxE_pv`](#--gxe-pv-arg)  
* [Visualisation](#--visualisation)
    * [`--kplots`](#--kplots-arg)
    * [`--distance`](#--distance-arg)
* [Additional Parameters](#Additional-Parameters)
    * [`--debug`](#--debug)
    * [`--version`](#--version)
    * [`--help`](#--help)
* [Software dependencies](#software-dependencies)
    * [`-profile`](#-profile)
    * [`-with-conda`](#-with-conda)
    * [`-with-docker`](#-with-docker)
    * [`-with-singularity`](#-with-singularity)
* [Other command line parameters](#other-command-line-parameters)
    * [`-work-dir`](#-work-dir)
    * [`-params-file`](#-params-file)
    * [`-config`](#-config)
    * [`-resume`](#-resume)
    * [`-name`](#-name)

## Workflow

![EpiDiverse/ewas Workflow](/docs/images/ewas_input.png)

## Running the pipeline
The main command for running the pipeline is as follows:

```bash
nextflow run epidiverse/ewas [OPTIONS]
```

Note that the pipeline will create files in your working directory:

```bash
work/           # Directory containing the nextflow working files
ewas/           # Finished results (configurable, see below)
.nextflow.log   # Log file from Nextflow
.nextflow/      # Nextflow cache and history information
```

## Inputs and Outputs

### `--input <ARG>` [REQUIRED]
Specify input path for the directory containing outputs from the WGBS pipeline. The pipeline searches for bedGraph files in '\*/bedGraph/{sample_name}_{context}.bedGraph' format, where sample names must correspond to the samplesheet and context can be either "CpG", "CHG", or "CHH".

### `--samples <ARG>` [REQUIRED]
Specify the path to the samplesheet file containing information regarding sample names and corresponding environment and covariate values. The file must contain at least three tab-separated columns: 1) sample names, 2) environment value, 3) covariate values, with further columns optional for additional covariates.

Example `samples.tsv` file:

```bash
#ID        env   cov1   cov2  . . . . .    covn
sample1	    34	  1      1
sample2	    42	  1      1  
sample3	    21	  1      2  
sample4	    56	  1      2         
sample5	    76	  2      3
sample6	    65	  2      3
sample7	    11	  2      4 
sample8	    22	  2      4
...
```
There can be multiple covariate columns but environment factor can only have one.

### `--DMPs <ARG>`
Specify path to the DMR pipeline output directory to run EWAS analyses in addition with methylated positions filtered by significant DMPs. The pipeline searches for bed files in '\*/{context}/metilene/\*/\*.bed' format where context can be either "CpG", "CHG", or "CHH".

### `--DMRs <ARG>`
Specify path to the DMR pipeline output directory to run EWAS analyses in addition with methylated positions filtered by significant DMRs. In addition, the pipeline will call the union of all significant regions and attempt to run EWAS with whole regions as markers. The pipeline searches for bed files in '*/{context}/metilene/*/*.bed' format where context can be either "CpG", "CHG", or "CHH".

### `--SNPs <ARG>`
ONLY SUITABLE FOR DIPLOID ORGANISMS. Specify path to the SNP pipeline output directory to enable EWAS analyses Gmodel and GxEmodel which attempt to create a genome-wide methQTL map. The pipeline searches for VCF files in '*/vcf/{sample_name}.{extension}' where sample names must correspond to the samplesheet and the extension can be any standard vcf extension readable by 'bcftools' and defined with --extension parameter. Alternatively, the path to a single multi-sample VCF file can be provided.

### `--extension <ARG>`
Specify the extension to use when searching for VCF files eg. \*.vcf \*.bcf or \*.vcf.gz [default: *.vcf.gz]

### `--output <ARG>`
A string that will be used as the name for the output results directory, which will be generated in the working directory. [default: ewas]


## Model Decision

### `--Emodel`
Run analysis with "E model". Disables other models unless they are also specified. If no individual model is specified then all that are possible with the provided inputs will run in parallel [default: off]

### `--Gmodel`
Run analysis with "G model". Disables other models unless they are also specified. If no individual model is specified then all that are possible with the provided inputs will run in parallel [default: off]

### `--GxE`
Run analysis with "GxE model". Disables other models unless they are also specified. If no individual model is specified then all that are possible with the provided inputs will run in parallel [default: off]

### `--noCpG`
Disables EWAS analysis in CpG context. Note: at least one methylation context is required for analysis. [default: off]

### `--noCHG`
Disables EWAS analysis in CHG context. Note: at least one methylation context is required for analysis. [default: off]

### `--noCHH`
Disables EWAS analysis in CHH context. Note: at least one methylation context is required for analysis. [default: off]


## Input Filtering

### `--coverage <ARG>`
Specify the minimum coverage threshold to filter individual methylated positions from the --input directory before running analyses [default: 0] 
            
### `--input_FDR <ARG>`
Specify the minimum FDR significance threshold to include DMPs and/or DMRs from the respective --DMPs and --DMRs directories [default: 0.05]          
                
### `--proportion <ARG>`
Minimum proportion of samples that must share a DMP and/or DMR for it to be considered in the analysis [default: 0.2]

### `--merge`
When running EWAS using the union set of DMRs as markers, specify to merge adjacent sub-regions into larger regions prior to methylation averaging and subsequent analysis [default: off]


## SNP Filtering 
NB: PROVIDING VARIANTS ONLY SUITABLE FOR DIPLOID ORGANISMS.

### `--max_missing <ARG>`
Variants that were successfully genotyped in given proportion of individuals. It can take values from 0 to 1, where 1 means no missing data allowed [default: 0.5]

### `--mac <ARG>`
Minor allele count [default: 3]

### `--minQ <ARG>`
Minimum quality score  [default: 30]


## Output Filtering
   
### `--Emodel_pv <ARG>`
Set the p-value to run "E model". Note: this filter is hardcoded as "1" for Q-Q plot generation and the user-given value is applied retroactively [default: 0.0001]

### `--Gmodel_pv <ARG>`
Set the p-value to run "G model". [default: 0.0001]

### `--GxE_pv <ARG>`
Set the p-value to run "GxE model". [default: 0.0001]

### `--output_FDR <ARG>`
Specify the maximum FDR threshold for filtering EWAS post-analysis [default: 0.05]


## Visualisation

### `--kplots <ARG>`
Specify the number of plots to generate for the top k significant results in "GxE model" [default: 10]

### `--distance <ARG>`
Specify the distance threshold to define cis and trans methQTLs in the dotplot generated for "G model" output [default: 5000]


## Additional Parameters

### `--debug`
Specify in order to prevent Nextflow from clearing the work dir cache following a successful pipeline completion. [default: off]

### `--version`
When called with `nextflow run epidiverse/ewas --version` this will display the pipeline version and quit.

### `--help`
When called with `nextflow run epidiverse/ewas --help` this will display the parameter options and quit.

## Software Dependencies

There are different ways to provide the required software dependencies for the pipeline. The recommended method is to use the Conda, Docker or Singularity profiles as provided by the pipeline. 

### `-profile`
Use this parameter to choose a preset configuration profile. See the [installation documentation](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/installation) for more information about profiles.

Profiles available with the pipeline are:

* `standard`
    * The default profile, used if `-profile` is not specified.
    * Uses sensible resource allocation for , runs using the `local` executor (native system calls) and expects all software to be installed and available on the `$PATH`.
    * This profile is mainly designed to be used as a starting point for other configurations and is inherited by most of the other profiles below.
* `conda`
    * Builds a conda environment from the environment.yml file provided by the pipeline
    * Requires conda to be installed on your system.
* `docker`
    * Launches a docker image pulled from epidiverse/dmr
    * Requires docker to be installed on your system. 
* `singularity`
    * Launches a singularity image pulled from epidiverse/dmr
    * Requires singularity to be installed on your system.
* `epi|diverse`
    * Designed to be used on the [EpiDiverse](http://epidiverse.eu/) clusters `epi` or `diverse`
    * Launches jobs using the `SLURM` executor.
    * Uses pre-built conda environments to provide all software requirements.
* `custom`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config for process resource allocation.

If you wish to provide your own package containers it is possible to do so by setting the `standard` or `custom` profile, and then providing your custom package with the command line flags below. These are not required with the the other profiles.

### `-with-conda <ARG>`
Flag to enable conda. You can provide either a pre-built environment or a *.yml file.

### `-with-docker <ARG>`
Flag to enable docker. The image will automatically be pulled from Dockerhub.

### `-with-singularity <ARG>`
Flag to enable use of singularity. The image will automatically be pulled from the internet. If running offline, follow the option with the full path to the image file.

## Other command line parameters

### `-work-dir <ARG>`
Specify the path to a custom work directory for the pipeline to run with (eg. on a scratch directory)

### `-params-file <ARG>`
Provide a file with specified parameters to avoid typing them out on the command line. This is useful for carrying out repeated analyses. A template params file [`assets/params.txt`](../assets/params.txt) has been made available in the pipeline repository.

### `-config <ARG>`
Provide a custom config file for adapting the pipeline to run on your own computing infrastructure. A template config file [`assets/custom.config`](../assets/custom.config) has been made available in the pipeline repository. This file can be used as a boilerplate for building your own custom config.

### `-resume [<ARG>]`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. Give a specific pipeline name as an argument to resume it, otherwise Nextflow will resume the most recent. NOTE: This will not work if the specified run finished successfully and the cache was automatically cleared. (see: [`--debug`](#--debug))

### `-name <ARG>`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
