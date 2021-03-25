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

1) Example of using Docker
```bash
nextflow run rikenbit/ramdaq -profile test,docker
```
2) Example of using Singularity
```bash
nextflow run rikenbit/ramdaq -profile test,singularity
```

iv. Start running your own analysis!

<!-- TODO: Update the default command above used to run the pipeline -->

iv-i. You can run ramdaq without donwloading reference annotation data.

```bash
nextflow run rikenbit/ramdaq -profile <docker/singularity> --reads '*_R{1,2}.fastq.gz' --genome GRCh38
```

iv-i. You can also run ramdaq by specifying local paths to reference annotation (See ['Using provided reference genome and annotations'](docs/local_annotation.md)).

```bash
nextflow run rikenbit/ramdaq -profile <docker/singularity> --reads '*_R{1,2}.fastq.gz' --genome GRCh38 --local_annot_dir <The directory path where the reference genome and annotations are placed>
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The ramdaq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. [Pipeline configuration](https://nf-co.re/usage/adding_own_config)
3. Running the pipeline
    * [Usage](docs/usage.md)
    * [Examples](docs/examples.md)
    * [Using test data](docs/test_data.md)
    * [Using bcl2fastq](https://github.com/rikenbit/ramdaq_bcl2fastq)
        * If you need to use BCL files produced by Illumina sequencing machines, execute [ramdaq_bcl2fastq](https://github.com/rikenbit/ramdaq_bcl2fastq).
        * [bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) is conversion software, which can be used to demultiplex data and convert BCL files to FASTQ file formats for downstream analysis.
        * Please see the README of [ramdaq_bcl2fastq](https://github.com/rikenbit/ramdaq_bcl2fastq) for details.
    * [Using provided reference genome and annotations](docs/local_annotation.md)
        * the current version supports human (GRCh38) and mouse (GRCm38).
    * [Using ramdaq on the NIG Supercomputer System](docs/nig_supercomputer_system.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits

ramdaq is written and maintained by [Mika Yoshimura](https://github.com/myoshimura080822) and [Haruka Ozaki](https://github.com/yuifu) in the collaboration of [Laboratory for Bioinformatics Research, RIKEN Center for Biosystems Dynamics Research](https://bit.riken.jp/) and [Bioinformatics Laboratory, Faculty of Medicine, University of Tsukuba](https://sites.google.com/view/ozakilab).

ramdaq was originally developed based on the [nf-core](https://nf-co.re/) template.

## Citation

[![DOI](https://zenodo.org/badge/269006630.svg)](https://zenodo.org/badge/latestdoi/269006630)
