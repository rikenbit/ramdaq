// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RSEQC  {
    tag "$name"
    label 'process_high'

    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true,
        saveAs: {filename ->
            if (filename.indexOf("readdist.txt") > 0)         "read_distribution/$filename"
            else if (filename.indexOf("inferexp.txt") > 0)    "infer_experiment/$filename"
            else if (filename.indexOf("inner_distance") > 0)  "inner_distance/$filename"
            else if (filename.indexOf("junction") > 0)        "junction_annotation/$filename"
            else if (filename.indexOf("splice_events") > 0)   "junction_annotation/$filename"
            else "$filename"
        }

    input:
    tuple val(name), file(bam), file(bai)
    file bed

    output:
    path "*.{txt,pdf,r,xls,log}", emit: rseqc_results
    path "${name}.readdist.txt", optional:true, emit: readdist_totalread

    script:
    def min_intron = params.min_intron > 0 ? "-m ${params.min_intron}" : ''
    if (params.single_end) {
        """
        read_distribution.py -i ${bam} -r $bed > ${bam.baseName}.readdist.txt
        infer_experiment.py -i ${bam} -r $bed > ${bam.baseName}.inferexp.txt
        junction_annotation.py -i ${bam} -o ${bam.baseName} -r $bed $min_intron 2> ${bam.baseName}.junction_annotation.log
        """
    } else {
        """
        read_distribution.py -i ${bam} -r $bed > ${bam.baseName}.readdist.txt
        infer_experiment.py -i ${bam} -r $bed > ${bam.baseName}.inferexp.txt
        junction_annotation.py -i ${bam} -o ${bam.baseName} -r $bed $min_intron 2> ${bam.baseName}.junction_annotation.log
        inner_distance.py -i ${bam} -o ${bam.baseName} -r $bed
        """
    }

}
