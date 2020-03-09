# EpiDiverse-EWAS Usage
This document describes the parameter options used by the pipeline.

* [Running the pipeline](#running-the-pipeline)
* [Inputs and outputs](#inputs-and-outputs)
    * [`--input`](#--input-arg-required)
    * [`--output`](#--output-arg)
    * [`--samples`](#--samples-arg-required)
    * [`--SNPs`](#--SNPs-arg-required-for-G-and-GXE-models)
    * [`--DMPs`](#--DMPs-arg)
    * [`--DMRs`](#--DMRs-arg)
    * [`--kplots`](#--kplots-arg-optional-for-GXE-model)
* [Inputs Filtering](#input-filtering)
    * [`--coverage`](#--coverage)
    * [`--input_FDR`](#--input_FDR)
    * [`--proportion`](#--proportion)
* [Model Decision](#Model-Decision)
    * [`--GEM_Emodel`](#--GEM_Emodel-arg-to-run-only-the-E-model)
    * [`--GEM_Gmodel`](#--GEM_Gmodel-arg-to-run-only-the-G-model)
    * [`--GEM_GXEmodel`](#--GEM_GXEmodel-arg-to-run-only-the-GXE-model) 
    * [`--GEM_Emodel --GEM_Gmodel`](#--GEM_Emodel--GEM_Gmodel-args-to-run-both-the-E-and-G-models)
    * [`--GEM_Emodel--GEM_GXEmodel`](#--GEM_Emodel--GEM_GXEmodel-args-to-run-both-the-E-and-GXE-models)
    * [`--GEM_Gmodel--GEM_GXEmodel`](#--GEM_Gmodel--GEM_GXEmodel-args-to-run-both-the-G-and-GXE-models)
    * [`--GEM_Emodel--GEM_Gmodel--GEM_GXEmodel`](#--GEM_Emodel--GEM_Gmodel--GEM_GXEmodel-args-to-run-all-models)
* [SNP Filtering](#--SNP-Filtering) 
    * [`--max_missing`](#--max_missing)
    * [`--mac`](#--mac)
    * [`--minQ`](#--minQ)
* [Output Filtering](#--Output-Filtering)
    * [`--output_FDR`](#--output_FDR)
* [Additional Parameters](#Additional-Parameters)
    * [`--noCpG`](#--noCpG)
    * [`--noCHG`](#--noCHG)
    * [`--noCHH`](#--noCHH)
    * [`--distance`](#--distance)
    * [`--extension`](#--extension)
    * [`--Emodel_pv`](#--extension)
    * [`--Gmodel_pv`](#--extension)
    * [`--GxEmodel_pv`](#--extension)    
    * [`--debug`](#--debug)
    * [`--version`](#--version)
    * [`--help`](#--help)
* [Software dependencies](#software-dependencies)
    * [`-profile`](#-profile)
    * [`-with-conda`](#-with-conda)
    * [`-with-docker`](#-with-docker)
    * [`-with-singularity`](#-with-singularity)    
## Workflow

![EpiDiverse/ewas Workflow](/docs/images/ewas_input.png)

## Inputs and Outputs

### `--input <ARG>` [REQUIRED]
Specify the path to the DIRECTORY containing each sample output from the WGBS pipeline to be taken forward for analysis. All the subdirectories must correspond to sample names in the provided samples.tsv file, and contain bedGraph directories with files in '*.bedGraph' format.

### `--output <ARG>`
A string that will be used as the name for the output results directory, which will be generated in the working directory. This directory will contain sub-directories for each set of reads analysed during the pipeline. [default: 'ewas']

### `--samples <ARG>` [REQUIRED]
Specify the path to the "samples.tsv" file containing information regarding sample names and corresponding groupings/replicates. The file must contain three tab-separated columns: 1) sample names, corresponding to subdirectories in the input directory, 1) sample names, 2) environment values and 3) covariate values

Example `*samples.tsv` file:

```bash
#sample_identifier          env   cov1   cov2  . . . . .    covn
FV_CZ_02_05_R1_GG0_M2_1	    34	  1      1
FV_CZ_02_05_R1_GG0_Y3_1	    42	  1      1  
FV_IT_01_06_R1_GG0_M2_1	    21	  1      2  
FV_IT_01_06_R1_GG0_Y3_1	    56	  1      2         
FV_IT_04_08_R1_GG0_M1_1	    76	  2      3
FV_IT_04_08_R1_GG0_Y2_1	    65	  2      3
FV_NO_01_06_R1_GG0_M2_1	    11	  2      4 
FV_NO_01_06_R1_GG0_Y3_1	    22	  2      4
...
```
Covariate data can contain multiple covariates but environment factor can only have one column.

### `--SNPs <ARG>` [required-for-G-and-GXE-models]
Specify the path to the DIRECTORY contaninig multiple .vcf (.bcf) outputs from the SNP pipeline. It is also possible to work with multi-sample merged vcf file. Pipeline can analyze both input types. The main aim is to convert genotypes to 1,2,3 matrix values for major allele homozygote (AA), heterozygote (AB) and minor allele homozygote (BB).

### `--DMPs <ARG>`
Specify the path to the DIRECTORY containing each sample output from the DMR pipeline to be taken forward for analysis with DMPs. All the subdirectories must correspond to sample names in the provided samples.tsv file, and contain bed file in '*.bed' format. It disables DMR run if it is used alone.

### `--DMRs <ARG>`
Specify the path to the DIRECTORY containing each sample output from the DMR pipeline to be taken forward for analysis with DMRs. All the subdirectories must correspond to sample names in the provided samples.tsv file, and contain bed file in '*.bed' format. It disables DMP run if it is used alone.

### `--kplots <ARG>` [optional-for-GXE-model]
Number of scatter plot to display ten top **significant** methylation corresponding to the environment in different genotype groups. AA, AB and BB are used for major allele homozygote, heterozygote and minor allele homozygote. Phenotypic values are visible on the x-axis, and methylation value in percentage on the y-axis. User can provide any number greater than zero. If any other methylations rather than significant ones want to be plotted, user should extract data manually and run it locally [default: 10]

## Input Filtering

### `--coverage`
Specify the minimum coverage threshold to filter methylated positions before running the EWAS analyses with all input types  [default: 0]

### `--input_FDR`
Specify the minimum coverage threshold to filter DMPs and/or DMR files  [default: 0.05]

### `--proportion`
Minimum proportion of samples that must share a DMP or DMR for it to be considered   [default: 0.20]


## Model Decision

### `--GEM_Emodel <ARG>` [to-run-only-the-E-model]

```bash
nextflow run epidiverse/ewas -profile <docker|singularity|conda> \
--input /path/to/reads/directory --samples /path/to/samples --GEM_Emodel
```

### `--GEM_Gmodel <ARG>` [to-run-only-the-G-model]


```bash
nextflow run epidiverse/ewas -profile <docker|singularity|conda> \
--input /path/to/reads/directory --samples /path/to/samples --SNPs /path/to/vcf/file(s) --GEM_Gmodel
```

### `--GEM_GXEmodel <ARG>` [to-run-only-the-GXE-model]

```bash
nextflow run epidiverse/ewas -profile <docker|singularity|conda> \
--input /path/to/reads/directory --samples /path/to/samples --SNPs /path/to/vcf/file(s) --GEM_Gmodel 
```

### `--GEM_Emodel--GEM_Gmodel <ARGS>` [to-run-both-the-E-and-G-models]

```bash
nextflow run epidiverse/ewas -profile <docker|singularity|conda> \
--input /path/to/reads/directory --samples /path/to/samples --SNPs /path/to/vcf/file(s) --GEM_Gmodel --GEM_Emodel
```

### `--GEM_Emodel--GEM_GXEmodel <ARGS>` [to-run-both-the-E-and-GXE-models]

```bash
nextflow run epidiverse/ewas -profile <docker|singularity|conda> \
--input /path/to/reads/directory --samples /path/to/samples --SNPs /path/to/vcf/file(s) --GEM_Emodel --GEM_GXEmodel
```

### `--GEM_Gmodel--GEM_GXEmodel <ARGS>` [to-run-both-the-G-and-GXE-models]

```bash
nextflow run epidiverse/ewas -profile <docker|singularity|conda> \
--input /path/to/reads/directory --samples /path/to/samples --SNPs /path/to/vcf/file(s) --GEM_Gmodel --GEM_GXEmodel
```

### `--GEM_Emodel--GEM_Gmodel--GEM_GXEmodel  <ARGS>` [--GEM_Emodel--GEM_Gmodel--GEM_GXEmodel-args-to-run-all-models] 

```bash
nextflow run epidiverse/ewas -profile <docker|singularity|conda> \
--input /path/to/reads/directory --samples /path/to/samples --SNPs /path/to/vcf/file(s) --GEM_Emodel --GEM_Gmodel --GEM_GXEmodel
```
or

```bash
nextflow run epidiverse/ewas -profile <docker|singularity|conda> \
--input /path/to/reads/directory --samples /path/to/samples --SNPs /path/to/vcf/file(s) 
```

Note that the pipeline will create files in your working directory:

```bash
work/           # Directory containing the nextflow working files
dmrs/           # Finished results (configurable, see below)
.nextflow.log   # Log file from Nextflow
.nextflow/      # Nextflow cache and history information
```


## SNP Filtering 
It is possible when user run the pipeline with --GEM_Gmodel parameter. Please Check this link for more details about fiiltering: http://www.ddocent.com/filtering/

### `--max_missing`
Variants that were successfully genotyped in given %[max missing] of individuals  [default: 0.5]. It can take values from 0 to 1, where 1 means no missing data allowed.

### `--mac`
Minor allele count [default: 3]

### `--minQ`
Minimum quality score  [default: 30]


## Output Filtering


### `--out_FDR`
Specify the maximum FDR threshold for filtering EWAS post-analysis.  [default: 0.05]

## Additional Parameters

### `--noCpG`
Disables EWAS analysis in CpG context. Note: at least one methylation context is required for analysis. [default: off]

### `--noCHG`
Disables EWAS analysis in CHG context. Note: at least one methylation context is required for analysis. [default: off]

### `--noCHH`
Disables EWAS analysis in CHH context. Note: at least one methylation context is required for analysis. [default: off]

### `--distance`
User can specify what defines cis and trans regions in the dotplot with Gmodel output. [default: 5000]

### `--extension`
User can use the pipeline with different input to the SNP pipeline can submit a directory with eg. .vcf .bcf or *.vcf.gz.

### `--Emodel_pv`
User can set the p-value to change blue and red line positions in Manhattan plot from Emodel run. [Default: 0.0001]
It is hardcoded as "1" for the Q-Q plot generation to use all p-values.

### `--Gmodel_pv`
User can set the p-value to run Gmodel. [Default: 0.0001]

### `--GxE_pv`
User can set the p-value to run GxEmodel. [Default: 0.0001]



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
