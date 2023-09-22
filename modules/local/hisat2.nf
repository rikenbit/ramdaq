// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process HISAT2  {
    tag "$name"
    label 'process_high'
    
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true,
        saveAs: { filename ->
                    filename.indexOf(".summary.txt") > 0 ? "summaries/$filename" : "$filename"
                }

    input:
    tuple val(name), file(reads)
    path hs2_indices
    path tools_dir

    output:
    tuple val(name), file("*.bam"), file("*.bai"), file("*.flagstat"), emit: hisat2_bam_all
    tuple val(name), file("${name}.bam"), file("${name}.bam.bai"), file("${name}.bam.flagstat"), optional:true, emit: hisat2_bam_qc
    tuple val(name), file("${name}.bam"), file("${name}.bam.bai"), file("${name}.bam.flagstat"), optional:true, emit: hisat2_bam_count
    path "*.summary.txt", emit: hisat2_summary

    script:
    def prefix = options.suffix ? "${name}${options.suffix}" : "${name}"

    def strandness = ''
    if (params.stranded == 'fr-firststrand') {
        strandness = params.single_end ? "--rna-strandness R" : "--rna-strandness RF"
    } else if (params.stranded == 'fr-secondstrand'){
        strandness = params.single_end ? "--rna-strandness F" : "--rna-strandness FR"
    }
    softclipping = params.softclipping ? '' : "--no-softclip"
    threads_num = params.hs_threads_num > 0 ? "-p ${params.hs_threads_num}" : ''
    index_base = hs2_indices[0].toString() - ~/.\d.ht2l?/

    if (params.single_end) {
        if (params.stranded && params.stranded != 'unstranded' && options.suffix != '.rrna') {
            """
            hisat2 $softclipping $threads_num -x $index_base -U $reads $strandness $options.args --summary-file ${prefix}.summary.txt \\
            | samtools view -bS - | samtools sort - -o ${prefix}.bam
            samtools index ${prefix}.bam
            samtools flagstat ${prefix}.bam > ${prefix}.bam.flagstat

            bamtools filter -in ${prefix}.bam -out ${prefix}.forward.bam -script ${tools_dir}/bamtools_f_SE.json
            samtools index ${prefix}.forward.bam
            samtools flagstat ${prefix}.forward.bam > ${prefix}.forward.bam.flagstat

            bamtools filter -in ${prefix}.bam -out ${prefix}.reverse.bam -script ${tools_dir}/bamtools_r_SE.json
            samtools index ${prefix}.reverse.bam
            samtools flagstat ${prefix}.reverse.bam > ${prefix}.reverse.bam.flagstat
            """
        } else {
            """
            hisat2 $softclipping $threads_num -x $index_base -U $reads $strandness $options.args --summary-file ${prefix}.summary.txt \\
            | samtools view -bS - | samtools sort - -o ${prefix}.bam
            samtools index ${prefix}.bam
            samtools flagstat ${prefix}.bam > ${prefix}.bam.flagstat
            """
        }

    } else {
        if (params.stranded && params.stranded != 'unstranded' && options.suffix != '.rrna') {
            """
            hisat2 $softclipping $threads_num -x $index_base -1 ${reads[0]} -2 ${reads[1]} $strandness $options.args --summary-file ${prefix}.summary.txt \\
            | samtools view -bS - | samtools sort - -o ${prefix}.bam
            samtools index ${prefix}.bam
            samtools flagstat ${prefix}.bam > ${prefix}.bam.flagstat

            bamtools filter -in ${prefix}.bam -out ${prefix}.forward.bam -script ${tools_dir}/bamtools_f_PE.json
            samtools index ${prefix}.forward.bam
            samtools flagstat ${prefix}.forward.bam > ${prefix}.forward.bam.flagstat

            bamtools filter -in ${prefix}.bam -out ${prefix}.reverse.bam -script ${tools_dir}/bamtools_r_PE.json
            samtools index ${prefix}.reverse.bam
            samtools flagstat ${prefix}.reverse.bam > ${prefix}.reverse.bam.flagstat

            samtools view -bS -f 0x40 ${prefix}.bam -o ${prefix}.R1.bam
            samtools index ${prefix}.R1.bam
            samtools flagstat ${prefix}.R1.bam > ${prefix}.R1.bam.flagstat

            samtools view -bS -f 0x80 ${prefix}.bam -o ${prefix}.R2.bam
            samtools index ${prefix}.R2.bam
            samtools flagstat ${prefix}.R2.bam > ${prefix}.R2.bam.flagstat
            """
        } else if (params.stranded == 'unstranded' && options.suffix != '.rrna'){
            """
            hisat2 $softclipping $threads_num -x $index_base -1 ${reads[0]} -2 ${reads[1]} $options.args --summary-file ${prefix}.summary.txt \\
            | samtools view -bS - | samtools sort - -o ${prefix}.bam
            samtools index ${prefix}.bam
            samtools flagstat ${prefix}.bam > ${prefix}.bam.flagstat

            samtools view -bS -f 0x40 ${prefix}.bam -o ${prefix}.R1.bam
            samtools index ${prefix}.R1.bam
            samtools flagstat ${prefix}.R1.bam > ${prefix}.R1.bam.flagstat

            samtools view -bS -f 0x80 ${prefix}.bam -o ${prefix}.R2.bam
            samtools index ${prefix}.R2.bam
            samtools flagstat ${prefix}.R2.bam > ${prefix}.R2.bam.flagstat
            """
        } else {
            """
            hisat2 $softclipping $threads_num -x $index_base -1 ${reads[0]} -2 ${reads[1]} $strandness $options.args --summary-file ${prefix}.summary.txt \\
            | samtools view -bS - | samtools sort - -o ${prefix}.bam
            samtools index ${prefix}.bam
            samtools flagstat ${prefix}.bam > ${prefix}.bam.flagstat
            """
        }
    }
}