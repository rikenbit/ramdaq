# CHANGELOG

## ramdaq: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## v1.1 (`-r 1.1`) [2021-03-30]

### `Added`

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
