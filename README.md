# ![ramdaq](docs/images/ramdaq_logo_vectorized.svg)

**This pipeline analyses data from full-length single-cell RNA sequencing (scRNA-seq) methods.**

[![GitHub Actions CI Status](https://github.com/rikenbit/ramdaq/workflows/CI/badge.svg)](https://github.com/rikenbit/ramdaq/actions)
[![GitHub Actions Linting Status](https://github.com/rikenbit/ramdaq/workflows/linting/badge.svg)](https://github.com/rikenbit/ramdaq/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

<!-- [![Docker](https://img.shields.io/docker/automated/nfcore/ramdaq.svg)](https://hub.docker.com/r/myoshimura080822/ramdaq) -->

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility (see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles)). Note that ramdaq does not support conda.

iii. Download the pipeline automatically and test it on a minimal dataset with a single command

```bash
nextflow run rikenbit/ramdaq -profile test,<docker/singularity>
```

> Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) and enable either `docker` or `singularity` to set the appropriate execution settings for your local computing environment.

iv. Start running your own analysis!

<!-- TODO: Update the default command above used to run the pipeline -->

```bash
nextflow run rikenbit/ramdaq -profile <docker/singularity> --reads '*_R{1,2}.fastq.gz' --genome GRCh37 --local_annot_dir <The directory path where the regerence genome and annotations are placed>
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The ramdaq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. [Pipeline configuration](https://nf-co.re/usage/adding_own_config)
3. Running the pipeline
    * [Usage](docs/usage.md)
    * [Using test data](docs/test_data.md)
    * [Using bcl2fastq](bcl2fastq/README.md)
        * [bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) is conversion software, which can be used to demultiplex data and convert BCL files to FASTQ file formats for downstream analysis.
    * [Using provided reference genome and annotations](docs/local_annotation.md)
        * the current version supports human (GRCh38) and mouse (GRCm38).
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

The ramdaq pipeline was originally developed based on the [nf-core](https://nf-co.re/) template.

<!-- TODO: Add a brief overview of what the pipeline does and how it works -->

## Credits

ramdaq was originally written by [Mika Yoshimura](https://github.com/myoshimura080822) and [Haruka Ozaki](https://github.com/yuifu) in the collaboration of [Laboratory for Bioinformatics Research, RIKEN Center for Biosystems Dynamics Research](https://bit.riken.jp/) and [Bioinformatics Laboratory, Faculty of Medicine, University of Tsukuba](https://sites.google.com/view/ozakilab).

## Citation

<!-- TODO: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  ramdaq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->You can cite the `nf-core` publication as follows:
