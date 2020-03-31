[<img width="200" align="right" src="docs/images/euflagbetter.jpg">](https://ec.europa.eu/programmes/horizon2020/en)
[<img width="200" align="right" src="docs/images/epidiverse-logo.jpg">](https://epidiverse.eu)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.09.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/epidiverse/ewas.svg)](https://hub.docker.com/r/epidiverse/ewas)

EpiDiverse-EWAS Pipeline
========================

**EpiDiverse/ewas** is a bioinformatics analysis pipeline for performing epigenome-wide association studies from methylated positions and/or regions, with optional analysis of methQTLs for diploid organisms from variant call data.

The workflow processes a population of sample bed files, usually derived from the [EpiDiverse-WGBS](https://github.org/epidiverse/wgbs) pipeline, and formats them with [bedtools](https://github.com/arq5x/bedtools2) for analysis with the R package [GEM](https://github.com/fastGEM/GEM). Output in the form of DMPs or DMRs from the [EpiDiverse-DMR](https://github.org/epidiverse/dmr) pipeline can also be given to pre-filter the number of positions and reduce multiple comparisons. In addition, the union set of DMRs can themselves be used as independent markers within EWAS in place of individual positions. Sample variants, usually derived from the [EpiDiverse-SNP](https://github.org/epidiverse/snp) pipeline, can additionally be provided to test the association of methylation-SNP pairs as methQTLs. The pipeline provides visualisation in the form of Manhattan plots for Emodel, sequence dotplots for Gmodel, genotype interaction plots for GxE model, and p-value QQ-plots + histograms for all models, each with [ggplot2](https://github.com/tidyverse/ggplot2).

> See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://www.nextflow.io/)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run epidiverse/ewas -profile test,<docker|singularity|conda>
```

iv. Start running your own analysis!

```bash
nextflow run epidiverse/ewas -profile <docker|singularity|conda> \
--input /path/to/wgbs/directory --samples /path/to/samples.tsv
```

> See the [usage documentation](docs/usage.md) for all of the available options when running the pipeline.

### Wiki Documentation

The EpiDiverse/template pipeline is part of the [EpiDiverse Toolkit](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/overview), a best practice suite of tools intended for the study of [Ecological Plant Epigenetics](https://app.gitbook.com/@epidiverse/s/project/). Links to general guidelines and pipeline-specific documentation can be found below:

1. [Installation](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/installation)
2. Pipeline configuration
    * [Local installation](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/installation#2-install-the-pipeline)
    * [Adding your own system config](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/installation#3-pipeline-configuration)
    * [EpiDiverse infrastructure](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/installation#appendices)
3. [Running the pipeline](docs/usage.md)
4. [Understanding the results](docs/output.md)
5. [Troubleshooting](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/troubleshooting)

### Credits

These scripts were originally written for use by the [EpiDiverse European Training Network](https://epidiverse.eu/), by Adam Nunn ([@bio15anu](https://github.com/bio15anu)) and Nilay Can ([@nilaycan](https://github.com/nilaycan)).

This project has received funding from the European Union’s Horizon 2020 research and innovation
programme under the Marie Skłodowska-Curie grant agreement No 764965

## Citation

If you use epidiverse/ewas for your analysis, please cite it using the following doi: <placeholder>
