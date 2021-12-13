// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALC_ERCC_CORRELATION {

    label 'process_low'
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true,
        saveAs: {filename ->
            filename.indexOf(".txt") > 0 ? "${options.publish_dir}/$filename" : "${options.args}/$filename"
        }
    
    input:
    val ercc_input_amount
    val num_of_files
    file(ercc_count)
    file(ercc_data)
    file(ercc_header)
    file(ercc_header_gstat)
    
    output:
    path "*.{txt,pdf}", emit: ercc_correlation_results
    path "barplot_*.csv", emit: ercc_correlation_barplot
    path "gstat_*.csv", emit: ercc_correlation_gstat
    
    //when:
    //num_of_files > 0

    script:
    """
    drawplot_ERCC_corr.r $ercc_count $ercc_data $ercc_input_amount
    cp barplot_ercc_counts_copynum_correlation.csv gstat_ercc_counts_copynum_correlation.csv
    cat $ercc_header barplot_ercc_counts_copynum_correlation.csv >> tmp_file
    mv tmp_file barplot_ercc_counts_copynum_correlation_mqc.csv
    cat $ercc_header_gstat gstat_ercc_counts_copynum_correlation.csv >> tmp_file
    mv tmp_file gstat_ercc_counts_copynum_correlation_mqc.csv
    """
}
