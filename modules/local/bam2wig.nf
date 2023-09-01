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
    def wigsum = params.wigsum && params.wigsum > 0 ? "--wigsum ${params.wigsum}" : ''
    """
    bam2wig.py -i ${bam} -s $chrsize ${wigsum} -u -o ${bam.baseName}
    """
}
