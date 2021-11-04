// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FASTQMCF  {
    tag "$name"

    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
 
    input:
    tuple val(name), file(reads)
    file adapter

    output:
    tuple val(name), file("*.fastq.gz"), emit: trimmed_reads

    script:
    maxReadLength = params.maxReadLength > 0 ? "-L ${params.maxReadLength}" : ''
    minReadLength = params.minReadLength > 0 ? "-l ${params.minReadLength}" : ''
    skew = params.skew > 0 ? "-k ${params.skew}" : ''
    quality = params.quality > 0 ? "-q ${params.quality}" : ''

    def prefix = options.suffix ? "${name}${options.suffix}" : "${name}"
    def prefix_1 = options.suffix ? "${name}_1${options.suffix}" : "${name}_1"
    def prefix_2 = options.suffix ? "${name}_2${options.suffix}" : "${name}_2"

    if (params.single_end) {
        """
        fastq-mcf $adapter $reads -o ${prefix}.fastq $maxReadLength $minReadLength $skew $quality ; gzip ${prefix}.fastq
        """
    } else {
        """
        fastq-mcf $adapter ${reads[0]} ${reads[1]} -o ${prefix_1}.fastq -o ${prefix_2}.fastq $maxReadLength $minReadLength $skew $quality ; gzip ${prefix_1}.fastq && gzip ${prefix_2}.fastq
        """
    }
}
