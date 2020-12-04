# ramdaq: Documentation

The ramdaq documentation is split into the following files:

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