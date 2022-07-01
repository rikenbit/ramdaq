// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALC_INITIALREADS_MAPRATE {

    label 'process_low'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
    
    input:
    file(seqcounts_merged)
    file(totalread_merged)
    file(initialread_maprate_header)
    file(initialread_maprate_header_gstat)
    
    output:
    path "*.{txt,pdf}", emit: initialread_maprate_results
    path "barplot_*.csv", emit: initialread_maprate_barplot
    path "gstat_initialread_maprate_*.csv", emit: initialread_maprate_gstat
    
    script:
    """
    drawplot_initialread_maprate_bar.r $seqcounts_merged $totalread_merged
    cp barplot_initialread_maprate.csv gstat_initialread_maprate.csv
    
    cat $initialread_maprate_header barplot_initialread_maprate.csv >> tmp_file
    mv tmp_file barplot_initialread_maprate_mqc.csv

    cat $initialread_maprate_header_gstat gstat_initialread_maprate.csv >> tmp_file
    mv tmp_file gstat_initialread_maprate_mqc.csv
    """
}
