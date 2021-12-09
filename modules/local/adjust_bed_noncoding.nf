// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ADJUST_BED_NONCODING {

    label 'process_low'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
    
    input:
    file(bed)
    
    output:
    path "adjusted.bed", emit: bed_adjusted
    
    script:
    """
    adjust_bed_noncoding.r $bed
    """
}
