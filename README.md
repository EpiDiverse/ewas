[<img width="200" align="right" src="docs/images/euflagbetter.jpg">](https://ec.europa.eu/programmes/horizon2020/en)
[<img width="200" align="right" src="docs/images/epidiverse-logo.jpg">](https://epidiverse.eu)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.09.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)

**EpiDiverse/ewas** is a bioinformatics analysis pipeline for aligning performing epigenome wide association studies (EWAS)for non-model plant species.

The workflow processes three input types, one is methylation calls from EpiDiverse WGBS pipeline ([WGBS] https://github.com/EpiDiverse/wgbs), others are Differentially Methylated Regions (DMRs) and Differentially Methylated Position (DMPs) from ([DMR] https://github.com/EpiDiverse/dmr) pipeline. EWAS analysis is performed by ([GEM] https://rdrr.io/bioc/GEM/man/GEM-package.html).
=======
[![Docker](https://img.shields.io/docker/automated/epidiverse/snp.svg)](https://hub.docker.com/r/epidiverse/snp)

EpiDiverse-template Pipeline
========================

**EpiDiverse/template** is a bioinformatics analysis pipeline for XXX with YYY data and this is a brief summary of the pipeline.

The workflow processes a collection of AAA files from the [Tool1](https://link/to/tool/1) pipeline using [Tool2](https://link/to/tool/2), then BBB and CCC. This paragraph gives a more involved description of the step-by-step processes that make up the pipeline. Finally the main output is described here which is performed using [Tool3](https://link/to/tool/3) and/or [Tool4](https://link/to/tool/4).

> See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://www.nextflow.io/)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run epidiverse/template -profile test,<docker|singularity|conda>
```

iv. Start running your own analysis!

```bash
nextflow run epidiverse/template -profile <docker|singularity|conda> --input /path/to/wgbs/directory
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

These scripts were originally written for use by the [EpiDiverse European Training Network](https://epidiverse.eu/), by Firstname Surname ([@username](https://github.com/link-to-user-profile)).

This project has received funding from the European Union’s Horizon 2020 research and innovation
programme under the Marie Skłodowska-Curie grant agreement No 764965

## Citation

If you use epidiverse/template for your analysis, please cite it using the following doi: <placeholder>
