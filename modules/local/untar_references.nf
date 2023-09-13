// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UNTAR_REFERENCES {

    label 'process_low'
    publishDir "${params.outdir}/${params.species}", mode: 'copy', overwrite: true
    
    input:
    path gz
    val common_dirname
    
    output:
    path "*", emit: untar_files
    
    script:
    def untar = gz.toString() - '.tar.gz'
    def mv_files = options.args ? "mv ${common_dirname}/* . && rmdir ${common_dirname}" : ''
    """
    tar -xvf $gz
    $mv_files
    """
}