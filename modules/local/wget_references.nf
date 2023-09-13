// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process WGET_REFERENCES {

    label 'process_low'
    //publishDir "${params.outdir}/${params.species}", mode: 'copy', overwrite: true
    
    input:
    val gz
    
    output:
    path "*.gz", emit: dl_file
    path "md5sum.txt", emit: md5sum
    
    script:
    fname = gz.toString()
    untar = gz.toString() - '.tar.gz'
    //tar -xvf $gz

    """
    wget ${params.dl_directory_url}/${params.species}/$fname &> /dev/null
    md5sum $fname | cut -c 1-32 > md5sum.txt
    """
}