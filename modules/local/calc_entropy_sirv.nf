// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALC_ENTROPY_SIRV {

    label 'process_low'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
    
    input:
    file(tpm_count)
    file(entropy_header)
    file(entropy_header_gstat)
    
    output:
    path "*.{txt,pdf}", emit: entropy_results
    path "barplot_*.csv", emit: entropy_barplot
    path "gstat_*.csv", emit: entropy_gstat

    script:
    def prefix = "entropy_${options.suffix}"
    """
    drawplot_sirv_entropy.r $tpm_count $prefix
    cp barplot_${prefix}.csv gstat_${prefix}.csv
    cat $entropy_header barplot_${prefix}.csv >> tmp_file
    mv tmp_file barplot_${prefix}_mqc.csv
    cat $entropy_header_gstat gstat_${prefix}.csv >> tmp_file
    mv tmp_file gstat_${prefix}_mqc.csv 
    """
}
