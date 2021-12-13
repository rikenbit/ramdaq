// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALC_SAMPLE_CORRELATION {

    label 'process_medium'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
    
    input:
    file(input_files)
    val num_of_bams
    file(mdsplot_header)
    file(heatmap_header)
    
    output:
    path "*.{txt,pdf,csv}", emit: sample_correlation
    
    when:
    num_of_bams > 2

    script:
    """
    edgeR_heatmap_MDS.r $input_files
    cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv >> tmp_file
    mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
    cat $heatmap_header log2CPM_sample_correlation_mqc.csv >> tmp_file
    mv tmp_file log2CPM_sample_correlation_mqc.csv
    """
}
