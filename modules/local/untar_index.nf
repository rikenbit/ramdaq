// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UNTAR_INDEX {

    label 'process_low'
    publishDir path: { params.saveReference ? "${params.outdir}/${options.publish_dir}" : params.outdir },
       saveAs: { params.saveReference ? it : null }, mode: 'copy'
    
    input:
    path gz
    
    output:
    path "$untar/*.ht2", emit: index_files
    
    script:
    untar = gz.toString() - '.tar.gz'
    """
    tar -xvf $gz
    """
}