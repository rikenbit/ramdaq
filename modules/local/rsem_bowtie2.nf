// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RSEM_BOWTIE2  {
    tag "$name"
    label 'process_high'
    
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true,
        saveAs: { filename ->
                    filename.indexOf(".log") > 0 ? "logs/$filename" : "$filename"
                }

    input:
    tuple val(name), file(bam), file(bai) //for input check
    tuple val(name), file(reads)
    path rsem_indices

    output:
    path "*.isoforms.results", emit: rsem_isoforms_to_merge
    path "*.genes.results", emit: rsem_genes_to_merge
    path "*.stat/*.cnt", emit: rsem_results_stat
    path "*.{results,log}"

    script:
    def prefix = options.suffix ? "${name}${options.suffix}" : "${name}"

    def strandness = ''
    if (params.stranded == 'fr-firststrand') {
        strandness = "--strandedness reverse"
    } else if (params.stranded == 'fr-secondstrand'){
        strandness = "--strandedness forward"
    }
    threads_num = params.rsem_threads_num > 0 ? "-p ${params.rsem_threads_num}" : ''
    index_base = rsem_indices[0].toString().split('\\.')[0]

    if (params.single_end) {
        if (params.stranded && params.stranded != 'unstranded') {
            """
            rsem-calculate-expression $threads_num $strandness $reads --bowtie2 --bowtie2-path /opt/conda/envs/ramdaq-1.0dev/bin/ $index_base ${prefix}
            samtools sort ${prefix}.transcript.bam -o ${prefix}.rsem.bam
            samtools index ${prefix}.rsem.bam
            samtools flagstat ${prefix}.rsem.bam > ${prefix}.rsem.bam.flagstat
            rm ${prefix}.transcript.bam
            """
        } else {
            """
            rsem-calculate-expression $threads_num $reads --bowtie2 --bowtie2-path /opt/conda/envs/ramdaq-1.0dev/bin/ $index_base ${prefix}
            samtools sort ${prefix}.transcript.bam -o ${prefix}.rsem.bam
            samtools index ${prefix}.rsem.bam
            samtools flagstat ${prefix}.rsem.bam > ${prefix}.rsem.bam.flagstat
            rm ${prefix}.transcript.bam
            """
        }
    } else {
        if (params.stranded && params.stranded != 'unstranded') {
            """
            rsem-calculate-expression $threads_num $strandness --paired-end ${reads[0]} ${reads[1]} --bowtie2 --bowtie2-path /opt/conda/envs/ramdaq-1.0dev/bin/ \\
            $index_base ${prefix}
            samtools sort ${prefix}.transcript.bam -o ${prefix}.rsem.bam
            samtools index ${prefix}.rsem.bam
            samtools flagstat ${prefix}.rsem.bam > ${prefix}.rsem.bam.flagstat
            rm ${prefix}.transcript.bam
            """
        } else {
            """
            rsem-calculate-expression $threads_num --paired-end ${reads[0]} ${reads[1]} --bowtie2 --bowtie2-path /opt/conda/envs/ramdaq-1.0dev/bin/ \\
            $index_base ${prefix}
            samtools sort ${prefix}.transcript.bam -o ${prefix}.rsem.bam
            samtools index ${prefix}.rsem.bam
            samtools flagstat ${prefix}.rsem.bam > ${prefix}.rsem.bam.flagstat
            rm ${prefix}.transcript.bam
            """
        }
    }
}