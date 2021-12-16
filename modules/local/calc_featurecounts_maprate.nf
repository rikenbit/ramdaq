// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALC_FEATURECOUNTS_MAPRATE {

    label 'process_low'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
    
    input:
    file(totalseq_merged)
    file(counts_merged)
    file(fcounts_maprate_header)
    file(fcounts_maprate_header_gstat)
    
    output:
    path "*.{txt,pdf}", emit: fcounts_maprate_results
    path "barplot_*.csv", emit: fcounts_maprate_barplot
    path "gstat_*.csv", emit: fcounts_maprate_gstat

    script:
    def prefix = "assignedrate_${options.suffix}"
    def is_pairedend = params.single_end ? "False" : "True"
    """
    drawplot_fcount_mappedrate_bar.r $totalseq_merged $counts_merged $is_pairedend ${options.suffix}
    cp barplot_${prefix}.csv gstat_${prefix}.csv
    cat $fcounts_maprate_header barplot_${prefix}.csv >> tmp_file
    mv tmp_file barplot_${prefix}_mqc.csv
    cat $fcounts_maprate_header_gstat gstat_${prefix}.csv >> tmp_file
    mv tmp_file gstat_${prefix}_mqc.csv
    """
}
