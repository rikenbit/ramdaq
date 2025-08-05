// Import generic module functions
include { initOptions } from './functions'

process MERGE_SUMMARYFILE {

    label 'process_low'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
    
    input:
    file(input_files)
    val options
    
    output:
    path "merged_*.txt", emit: merged_summary
    
    script:
    def prefix = "merged_${options.suffix}"
    command = input_files.collect{filename ->
        "awk '{if (FNR==1){print FILENAME, FNR, NR, \$0}}' ${filename} | sed '$options.args' | cut $options.args2 --delim=\" \" >> ${prefix}.txt"}.join(" && ")
    """
    $command
    """
}
