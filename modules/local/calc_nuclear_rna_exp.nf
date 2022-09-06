// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALC_NUCLEAR_RNA_EXP {

    label 'process_low'
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true,
        saveAs: {filename ->
            filename.indexOf(".txt") > 0 ? "${options.publish_dir}/$filename" : "${options.args}/$filename"
        }
    
    input:
    file(raw_count)
    file(tpm_count)
    file(plot_header)
    
    output:
    path '*_nuclearRNA_log.txt', optional:true, emit: nuclearrna_tpm_merged
    path "barplot_*.csv", emit: nuclear_rna_barplot
    
    script:
    def nuclearRNA_qc = params.nuclearRNA_qc ? params.nuclearRNA_qc : 'FALSE'
    
    """
    drawplot_nuclear_rna_exp.r $raw_count $tpm_count $nuclearRNA_qc
    cat $plot_header barplot_nuclear_rna_exp.csv >> tmp_file
    mv tmp_file barplot_nuclear_rna_exp_mqc.csv
    """
}
