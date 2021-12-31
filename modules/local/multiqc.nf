// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MULTIQC {

    label 'process_medium'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/multiqc:1.10.1--pyhdfd78af_1"

    input:
    path multiqc_config
    path multiqc_custom_config
    path ('fastqc/*')
    path ('fastqc/*')
    path ('hisat2_genome/*')
    path ('hisat2_rrna/*')
    path ('rseqc/*')
    path ('featureCounts/biotype_counts/*')
    path ('rsem_bowtie2_allgenes/*')
    path ('plot_sample_correlation/*')
    path ('plot_ercc_correlation/*')
    path ('plot_ercc_correlation/*')
    path ('featurecounts_all_gtf/*')
    path ('featurecounts_mt_gtf/*')
    path ('featurecounts_histone_gtf/*')
    path ('plot_assignedgenome/*')
    path ('plot_assignedgenome/*')
    path ('plot_fcounts_maprate_allgene/*')
    path ('plot_fcounts_maprate_allgene/*')
    path ('plot_fcounts_maprate_mt/*')
    path ('plot_fcounts_maprate_mt/*')
    path ('plot_fcounts_maprate_histone/*')
    path ('plot_fcounts_maprate_histone/*')
    path ('plot_detectedgenes_dr/*')
    path ('plot_detectedgenes_dr/*')
    path ('plot_detectedgenes_dr/*')
    path ('plots_from_tpmcounts_rsem/*')
    path ('plots_from_tpmcounts_rsem/*')
    path ('plots_from_tpmcounts_rsem/*')
    path ('plots_from_tpmcounts_rsem/*')
    path ('plots_entropy_sirv/*')
    path ('plots_entropy_sirv/*')
    path ('software_versions/*')
    path workflow_summary
    
    output:
    path "*multiqc_report.html", emit: multiqc_report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    
    script:
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc -f $custom_config .
    """
}