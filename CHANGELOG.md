# CHANGELOG

## ramdaq: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## v1.9 (`-r 1.9`) [2022-07-20]

### `v1.9 Added`

- Add barplot alligned genome rates of initial fastq reads

### `v1.9 Fixed`

- Change setting testcase resource limits
- Fixed the unintended samples omit problem in low mapped reads bam QC processes
- Remove duplicate column from General statistics
- Changed bamtools filter option for create strand-specific bams

## v1.8 (`-r 1.8`) [2022-06-21]

### `v1.8 Fixed`

- Fixed a problem with Fastq-mcf incorrect arguments that prevented fastq trimming from working properly

## v1.7 (`-r 1.7`) [2022-05-19]

### `v1.7 Added`

- Support for genocode v39/vM28, v40/vM29 annotations

### `v1.7 Fixed`

- Fixed an error in the RSEM-SIRV process when execute PE unstranded mode
- Extended "max time" of process-specific resource requirements
- Added hidden option to changed SIRV coverage process to off by default
- Added hidden option to avoid error when only short introns are detected in RSeQC

## v1.6 (`-r 1.6`) [2022-02-22]

- Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)

### `v1.6 Added`

- Add troubleshooting.md

### `v1.6 Fixed`

- Fixed an error when running only one sample

## v1.5 (`-r 1.5`) [2021-12-20]

### `v1.5 Added`

- Add surge.sh publish action to all CI test cases
- Support for genocode v38/vM27 annotations
- Update the docker image (some technical changes)

## v1.4.1 (`-r 1.4.1`) [2021-09-09]

### `v1.4.1 Fixed`

- Fixed the problem that dimensional reduction plots are not displayed in the MulitiQC report

## v1.4 (`-r 1.4`) [2021-08-23]

### `v1.4 Added`

- Added new rRNA quantification process and multiQC plot. More accurate % rRNA is achieved by mapping reads on rRNA sequences using HISAT2.
- Gzip-compressed output of `merge_readcoverage_sirv`

## v1.3 (`-r 1.3`) [2021-07-29]

### `v1.3 Added`

- Added SIRV entropy bar plot as new SIRV QC plot
- Added columns to the General Statistics table for all bar plots
- Added the summary plot of featurecounts in MultiQC

## v1.2.2 (`-r 1.2.2`) [2021-06-22]

### `v1.2.2 Fixed`

- Fixed the problem where hyphens in sample names were converted in MultiQC

## v1.2.1 (`-r 1.2.1`) [2021-06-14]

### `v1.2.1 Added`

- Added `--entire_max_cpus` and `--entire_max_memory` options to set upper limit of CPU and memory usage of the entire workflow.

### `v1.2.1 Fixed`

- Fixed the problem of high load average for some alignment steps
- Fixed a typo in the name of the log file that is copied to the output directory

## v1.2 (`-r 1.2`) [2021-05-26]

### `v1.2 Added`

- Support for outputting gene-level quantification results of RSEM
- Added RSEM results of the number of detected genes and transcripts in MultiQC report
- Changed the unit of the number of reads from M Seqs to K Seqs in MultiQC report
- New bigWig file output from R1, R2, forward, and reverse BAM files
- Save a copy of `.nextflow.log` in the output directory

### `v1.2 Fixed`

- Fixed the caluclation of read assgined rate of all genes, mitochondrial, rRNA, and histone by featureCounts
- Fixed rRNA annotation for GRCm39
- Fixed TPM calculation for ERCC
- Corrected URL of remote annotation files in `conf/remote_annotation.config`
- Improved human rRNA annotation GTF

## v1.1 (`-r 1.1`) [2021-03-30]

### `v1.1 Added`

- The new CI-test pattern for github actions
- A new option to allow the number of ERCC copies to be specified arbitrarily
- The process of transcripts quantification using RSEM-Bowtie2
- A new annotation files (gencode mouse vM26, human v37)
- The ability to run a quantitative process for each transcript if samples contain Spike-In RNA Variant (SIRV) Control
- A new bamQC plot "Junction Annotation" from RSeQC
- The process of adjusting BED in the workflow to account for noncoding region assignments to correct
- A new config "remote_annotation.config" obtaining annotations from [public data](https://bioinformatics.riken.jp/ramdaq)

### `Fixed`

- Markdown linting error in github actions
- The problem of "_1" and "_2" in the file name being erased in MultiQC report
- The problem of MultiQC causing memory error
- An error that occurred frequently at the end of a long pipeline
- An error that occurred when drawing a histone assigned plot at only one sample.
- Hide the aggregate results of rRNA in "General Stats" of MultiQC report
- Corrected the notation "assigned to gene" to "assigned to genome"
- Output the plot even if the data are all 0

### `Dependencies`

### `Deprecated`

- igenomes.config
