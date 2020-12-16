# ramdaq: Documentation

The ramdaq documentation is split into the following files:

1. [Installation](https://nf-co.re/usage/installation)
2. [Pipeline configuration](https://nf-co.re/usage/adding_own_config)
3. Running the pipeline
    * [Usage](usage.md)
    * [Using test data](test_data.md)
    * [Using bcl2fastq](https://github.com/rikenbit/ramdaq_bcl2fastq)
        * If you need to use BCL files produced by Illumina sequencing machines, execute [ramdaq_bcl2fastq](https://github.com/rikenbit/ramdaq_bcl2fastq).
        * [bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) is conversion software, which can be used to demultiplex data and convert BCL files to FASTQ file formats for downstream analysis.
        * Please see the README of [ramdaq_bcl2fastq](https://github.com/rikenbit/ramdaq_bcl2fastq) for details.
    * [Using provided reference genome and annotations](local_annotation.md)
        * the current version supports human (GRCh38) and mouse (GRCm38).
4. [Output and how to interpret the results](output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)
