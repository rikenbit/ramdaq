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
    path "merged_featureCounts_${options.suffix}.txt", emit: merged_counts
    path '*_TPM.txt', optional:true, emit: counts_tpm_merged
    path '*_ERCC_TPM_log.txt', optional:true, emit: ercc_tpm_merged
    
    script:
    def prefix = "merged_featureCounts_${options.suffix}"
    def get_gene_ids = "<(tail $options.args ${input_files[0]} | cut $options.args2 )"
    def get_counts = input_files.collect{filename ->
        "<(tail $options.args ${filename} | sed '1s/.bam//' | cut $options.args3)"}.join(" ")

    if (options.suffix != 'allgene'){
        """
        paste $get_gene_ids $get_counts > ${prefix}.txt
        """
    } else {
        """
        paste $get_gene_ids $get_counts > ${prefix}.txt
        calc_TPMCounts.r ${prefix}.txt
        """
    }

}
