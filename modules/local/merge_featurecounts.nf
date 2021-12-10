// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_FEATURECOUNTS {

    label 'process_low'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
    
    input:
    file(input_files)
    
    output:
    path "merged_*.txt", emit: merged_counts
    
    script:
    def prefix = "merged_featureCounts_${options.suffix}"
    def get_gene_ids = "<(tail $options.args ${input_files[0]} | cut $options.args2 )"
    def get_counts = input_files.collect{filename ->
        "<(tail $options.args ${filename} | sed 's:.bam::' | cut $options.args3)"}.join(" ")

    """
    paste $get_gene_ids $get_counts > ${prefix}.txt
    """
}
