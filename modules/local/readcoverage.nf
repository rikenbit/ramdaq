// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process READCOVERAGE  {
    tag "$name"
    label 'process_high'
    container "yuifu/readcoverage.jl:0.1.2-workaround"

    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true,
            saveAs: {filename ->
            if (filename.indexOf("geneBodyCoverage.txt") > 0) "genebody_coverage/$filename"
            else "$filename"
            }
 
    input:
    tuple val(name), file(bam), file(bai)
    file bed

    output:
    file "*.txt"

    script:
    """
    julia /opt/run.jl relcov ${bam} $bed ${bam.baseName}
    """
}
