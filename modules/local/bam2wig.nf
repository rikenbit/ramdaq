// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BAM2WIG  {
    label 'process_medium'
    tag "$name"

    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
 
    input:
    tuple val(name), file(bam), file(bai)
    file chrsize

    output:
    file "*.bw"
    file "*.wig"

    script:
    """
    bam2wig.py -i ${bam} -s $chrsize -u -o ${bam.baseName}
    """
}
