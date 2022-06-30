// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALC_ASSIGNEDGENOME_RATE {

    label 'process_low'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
    
    input:
    file(totalseq_merged)
    file(totalread_merged)
    file(assignedgenome_header)
    
    output:
    path "*.{txt,pdf}", emit: assignedgenome_rate_results
    path "barplot_*.csv", emit: assignedgenome_rate_barplot
    
    script:
    def is_pairedend = params.single_end ? "False" : "True"
    """
    drawplot_assignedgenomerate_bar.r $totalseq_merged $totalread_merged $is_pairedend
    cp barplot_assignedgenome_rate.csv gstat_assignedgenome_rate.csv
    cat $assignedgenome_header barplot_assignedgenome_rate.csv >> tmp_file
    mv tmp_file barplot_assignedgenome_rate_mqc.csv
    """
}
