#!/usr/bin/env nextflow

out_dir = params.out_dir
run_basedir = params.run_basedir

process bcl2fastq  {
    publishDir "${out_dir}/fastq_files"

    container "bcl2fastq2:1.0"

    input:
    val out_dir
    path run_basedir

    script:
    def output_d = "${out_dir}/fastq_files" 

    """
    bcl2fastq --no-lane-splitting --runfolder-dir $run_basedir --interop-dir $run_basedir/InterOp --input-dir $run_basedir/Data/Intensities/BaseCalls --sample-sheet $run_basedir/SampleSheet.csv --output-dir $PWD/$output_d --stats-dir $PWD/$output_d/Stats --reports-dir $PWD/$output_d/Reports && rm -rf $PWD/$output_d/Undetermined*
    """
}
