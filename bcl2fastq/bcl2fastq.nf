#!/usr/bin/env nextflow

proj_id = params.proj_id

Channel.from(params.run_basedirs)
    .into{run_basedirs; run_basedir_print;}

run_basedir_print.println()

process bcl2fastq  {
    tag { "${run_id}"}
    publishDir "output_${proj_id}/${run_id}/01_fastq_files"

    container "docker.io/myoshimura080822/bcl2fastq2:2.0"

    input:
    val proj_id
    set run_id, path(run_basedir) from run_basedirs

    script:
    def output_dir = "output_${proj_id}/${run_id}/01_fastq_files" 

    """
    bcl2fastq --no-lane-splitting --runfolder-dir $run_basedir --interop-dir $run_basedir/InterOp --input-dir $run_basedir/Data/Intensities/BaseCalls --sample-sheet $run_basedir/SampleSheet.csv --output-dir $PWD/$output_dir --stats-dir $PWD/$output_dir/Stats --reports-dir $PWD/$output_dir/Reports && rm -rf $PWD/$output_dir/Undetermined*
    """
}
