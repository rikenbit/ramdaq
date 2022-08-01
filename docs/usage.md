# ramdaq: Usage

## Table of contents

- [ramdaq: Usage](#ramdaq-usage)
  - [Table of contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Running the pipeline](#running-the-pipeline)
    - [Updating the pipeline](#updating-the-pipeline)
  - [Main arguments](#main-arguments)
    - [`-profile`](#-profile)
    - [`-c`](#-c)
    - [`--reads`](#--reads)
    - [`--single_end`](#--single_end)
    - [`--stranded`](#--stranded)
  - [Reference genomes and annotations](#reference-genomes-and-annotations)
    - [`--genome`](#--genome)
    - [`--saveReference`](#--savereference)
    - [`--local_annot_dir`](#--local_annot_dir)
    - [`--spike_in_ercc`](#--spike_in_ercc)
    - [`--spike_in_sirv`](#--spike_in_sirv)
  - [Other command line parameters](#other-command-line-parameters)
    - [`--outdir`](#--outdir)
    - [`-name`](#-name)
    - [`-resume`](#-resume)
    - [`--max_memory`](#--max_memory)
    - [`--max_time`](#--max_time)
    - [`--max_cpus`](#--max_cpus)
    - [`--monochrome_logs`](#--monochrome_logs)
  - [Parameters for each tools](#parameters-for-each-tools)
    - [Fastqmcf](#fastqmcf)
    - [Hisat2](#hisat2)
    - [featureCounts](#featurecounts)
  - [Results report options](#results-report-options)
    - [`--sampleLevel`](#--samplelevel)

<!--
    - [`--email`](#--email)
    - [`--email_on_fail`](#--email_on_fail)
    - [`--max_multiqc_email_size`](#--max_multiqc_email_size)
    - [`--plaintext_email`](#--plaintext_email)
    - [`--multiqc_config`](#--multiqc_config)
-->

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

<!-- TODO : Document required command line parameters to run the pipeline-->

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run rikenbit/ramdaq --reads '*_R{1,2}.fastq.gz' -profile docker  --genome GRCm39_vM26 --local_annot_dir <annotation directory path>
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull rikenbit/ramdaq
```

<!-- -
### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [ramdaq releases page](https://github.com/rikenbit/ramdaq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.
 -->

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](http://docker.com/)
  - Pulls software from dockerhub: [`myoshimura080822/ramdaq`](http://hub.docker.com/r/myoshimura080822/ramdaq/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  - Pulls software from DockerHub: [`myoshimura080822/ramdaq`](http://hub.docker.com/r/myoshimura080822/ramdaq/)
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_R{1,2}.fastq.gz'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.
4. The FASTQ files have to be 'gzipped' (e.g., `.fastq.gz`, `.fq.gz`). (Otherwise, You would get an error).

If left unspecified, a default pattern is used: `fastq_files/*{1,2}.fastq.gz`

### `--single_end`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--single_end --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

### `--stranded`

By default, the pipeline expects unstranded data. If you use the library preparation kit for NGS with orientation information, you need to specify the orientation of fastq in the stranded option. For example:

```bash
--stranded fr-firststrand
```

- unstranded : default (If not specified setting)
- fr-firststrand : First read corresponds to the reverse complemented counterpart of a transcript
- fr-secondstrand : First read corresponds to a transcript

## Reference genomes and annotations

### `--genome`

Specifies reference genome version. The available versions and specified options are as follows.

- Human
  - GRCh38
    - gencode.v37.primary_assembly.annotation: `--genome GRCh38_v37`
    - gencode.v35.primary_assembly.annotation: `--genome GRCh38_v35`
- Mouse
  - GRCm39
    - gencode.vM26.primary_assembly.annotation: `--genome GRCm39_vM26`
  - GRCm38
    - gencode.vM25.primary_assembly.annotation: `--genome GRCm38_vM25`

> For the implementation of the workflow, you need to download some mapping index or annotation files. Check the [this document](docs/local_annotation.md) for details.

Note that you can use the same configuration setup to save sets of reference files for your own use. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

<!-- TODO: Update reference genome example according to what is needed -->

```nextflow
params {
  genomes {
    'GRCh38' {
      gtf   = '<path to the gtf file>' // Used if no star index given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### `--saveReference`

Save the downloded reference files to the results directory

### `--local_annot_dir`

Specifes the location of local path to a directiory containing annotation files. See [local_annotation](local_annotation.md)

### `--spike_in_ercc`

Dilution rate of the ERCC Spike-In Control Mix 1 (default: false). Use when the samples contain the [ERCC Spike-In Control Mix 1](https://www.thermofisher.com/order/catalog/product/4456740#/4456740). The value is used to calculate the copy number of ERCC. If the value is not specified, '2e-7' is used as dilution rate. eg. `--spike_in_ercc '2e-7'`

When this option is specified, ramdaq outputs the between-sample correlation of expression levels of ERCC spike-ins in the MultiQC report.

### `--spike_in_sirv`

Dilution rate of the SIRV-Set 4 (default: false). Use when the samples contain the [SIRV-Set 4 (Lexogen's Spike-In RNA Variant Control)](https://www.lexogen.com/sirvs/). The value is used to calculate the copy number of ERCC in the SIRV-Set 4. eg. `--spike_in_sirv '4e-6'`

When this option is specified, ramdaq outputs the between-sample correlation of expression levels of ERCC spike-ins in the MultiQC report. Besides, the read coverages on the SIRV genes are saved in the output directory for downstream use.

## Computational performance

### Resource allocation for the entire workflow

#### `--entire_max_cpus`

Maximum number of CPUs to use for each step of the pipeline.
Should be a string in the format integer-unit. eg. `--entire_max_cpus 16`.

#### `--entire_max_memory`

Memory limit for each step of the pipeline.
Should be a string in the format integer-unit. eg. `--entire_max_memory '16.GB'`.

### Resource allocation for each process

#### `--max_cpus`

Use to set a top limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

#### `--max_memory`

Use to set a top limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

#### `--max_time`

Use to set a top limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

## Other command line parameters

<!-- TODO: Describe any other command line flags here -->

### `--outdir`

The output directory where the results will be saved.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

<!--
### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.

-->

## Parameters for each tools

### Fastqmcf

> [Fastqmcf](https://expressionanalysis.github.io/ea-utils/) : Scans a sequence file for adapters, and, based on a log-scaled threshold, determines a set of clipping parameters and performs clipping. Also does skewing detection and quality filtering.

- fastq-mcf attempt to detect and remove adapter and primer sequences.
  - ramdaq currently uses tge [fasta file](https://bioinformatics.riken.jp/ramdaq/ramdaq_annotation/mouse/all_sequencing_WTA_adopters.fa) that describe multiple adapter sequences.

- `--maxReadLength [N]`
  - Maximum remaining sequence length (Default: 75)
- `--minReadLength [N]`
  - Minimum remaining sequence length (Default: 36)
- `--skew [N]`
  - Skew percentage-less-than causing cycle removal (Default: 4)
- `--quality [N]`
  - Quality threshold causing base removal (Default: 30)

### Hisat2

> [HISAT2](http://daehwankimlab.github.io/hisat2/) : A fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA) to a population of human genomes as well as to a single reference genome.

- `--softclipping`
  - HISAT2 allow soft-clip reads near their 5' and 3' ends (Default: disallow)
- `--hs_threads_num [N]`
  - HISAT2 to launch a specified number of parallel search threads (Default: 1)

### featureCounts

> [featureCounts](http://subread.sourceforge.net/) : a software program developed for counting reads to genomic features such as genes, exons, promoters and genomic bins.

- `--extra_attributes`
  - Define which extra parameters should also be included in featureCounts (Default: 'gene_name')
- `--group_features`
  - Define the attribute type used to group features (Default: 'gene_id')
- `--count_type`
  - Define the type used to assign reads (Default: 'exon')
- `--allow_multimap`
  - Multi-mapping reads/fragments will be counted (Default: true)
- `--allow_overlap`
  - Reads will be allowed to be assigned to more than one matched meta-feature (Default: true)
- `--count_fractionally`
  - Assign fractional counts to features  (Default: true / This option must be used together with ‘--allow_multimap’ or ‘--allow_overlap’ or both)
- `--fc_threads_num [N]`
  - Number of the threads (Default: 1)
- `--group_features_type`
  - Define the type attribute used to group features based on the group attribute (default: 'gene_type')

## Results report options

### `--sampleLevel`

Used to turn off the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples.
