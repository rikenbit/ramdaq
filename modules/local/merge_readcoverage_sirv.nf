// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_READCOVERAGE_SIRV {

    label 'process_high'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
    
    input:
    file(input_files)
    
    output:
    path "*.gz", emit: coverage_sirv_merged
    
    script:
    def prefix = "merged_${options.suffix}.tsv"
    def command = input_files.collect{filename ->
        "awk -v filebasename=${filename.name} 'BEGIN{OFS=\"\t\"; bam=filebasename; sirv=filebasename; sub(\"\\..+?.readCoverage.txt\", \"\", bam); sub(\".readCoverage.txt\", \"\", sirv); sub(\".+\\.\", \"\", sirv);}{print sirv, \$1, \$2, bam}' ${filename} >> ${prefix}"}.join(" && ")

    """
    echo -e "SIRV\tposition\tcoverage\tbam" > ${prefix}
    $command
    gzip ${prefix}
    """
}
