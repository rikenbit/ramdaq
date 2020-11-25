# Nextflow pipeline for bcl2fastq

## Converting FASTQ files from a BCL file
First, pull a Docker container for bcl2fastq with the following command.
```
cd rikenbit/ramdaq.nf/bcl2fastq
docker build -t bcl2fastq2:1.0 .
```

Then, run test.
```
$ docker run --rm bcl2fastq2:1.0 bcl2fastq --help
BCL to FASTQ file converter
bcl2fastq v2.17.1.14
Copyright (c) 2007-2015 Illumina, Inc.
...
```
Then, follow the below instructions.

```
# Make a working directory with your favorite name and move into the directory
mkdir XXXXXX;cd XXXXXX

# Copy a config file (The example here is in case of SE)
cp rikenbit/ramdaq.nf/bcl2fastq/bcl2fastq.config .

# Rewrite config (with vi, vim, emacs, or your favorite editors)
# Setting requires a BaseCalls directory path(full) containing the binary base call files (BCL files) and output directory name.
vi bcl2fastq.config

# Run bcl2fastq via nextflow
nextflow run rikenbit/ramdaq.nf/bcl2fastq/bcl2fastq.nf -c bcl2fastq.config
```