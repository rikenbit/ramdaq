// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FASTQC {
    label 'process_low'
    tag "$name"
    
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy',
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    tuple val(name), file(reads)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    def prefix = options.suffix ? "${name}${options.suffix}" : "${name}"
    def prefix_1 = options.suffix ? "${name}_1${options.suffix}" : "${name}_1"
    def prefix_2 = options.suffix ? "${name}_2${options.suffix}" : "${name}_2"

    if (params.single_end) {
        newfastq = (reads.getName() =~ /\.gz$/) ? "${prefix}.fastq.gz" : "${prefix}.fastq"
        """
        if [ ! -f $newfastq ]; then
            ln -s $reads $newfastq
        fi
        fastqc $options.args --threads $task.cpus $newfastq
        """
    } else {
        newfastq1 = (reads[0].getName() =~ /\.gz$/) ? "${prefix_1}.fastq.gz" : "${prefix_1}}.fastq"
        newfastq2 = (reads[1].getName() =~ /\.gz$/) ? "${prefix_2}.fastq.gz" : "${prefix_2}.fastq"
        """
        if [ ! -f $newfastq1 ]; then
            ln -s ${reads[0]} $newfastq1
            ln -s ${reads[1]} $newfastq2
        fi
        fastqc $options.args --threads $task.cpus $newfastq1
        fastqc $options.args --threads $task.cpus $newfastq2
        """
    }
}
