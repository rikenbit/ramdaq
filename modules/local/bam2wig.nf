// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BAM2WIG  {
    //label 'process_low'
    tag "$name"

    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
 
    input:
    tuple val(name), file(bam), file(bai)
    file chrsize

    output:
    file "*.bw"

    script:
    """
    bam2wig.py -i ${bam} -s $chrsize --wigsum 100000000 -u -o ${bam.baseName}
    rm ${bam.baseName}.wig
    """
}
