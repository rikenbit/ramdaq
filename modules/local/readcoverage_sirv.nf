// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process READCOVERAGE_SIRV  {
    tag "$name"
    label 'process_high'
    container "yuifu/readcoverage.jl:0.1.2-workaround"

    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
 
    when:
    params.sirv_coverage

    input:
    tuple val(name), file(bam), file(bai)

    output:
    path "*.txt", emit: readcov_sirv_results

    script:
    def leftpos = [1001, 14644, 22555, 34498, 51620, 67226, 81063, 230020, 235017, 240016, 245017, 252017, 259017, 266017, 275018, 284018, 293018, 303959, 314960, 325930, 338959, 351958]
    def rightpos = [11643, 19554, 31497, 48619, 64225, 78062, 228019, 234016, 239015, 244016, 251016, 258016, 265016, 274017, 283017, 292017, 302958, 313959, 324929, 337958, 350957, 363957]
    def sirvname = ["SIRV1", "SIRV2", "SIRV3", "SIRV4", "SIRV5", "SIRV6", "SIRV7", "SIRV4001", "SIRV4002", "SIRV4003", "SIRV6001", "SIRV6002", "SIRV6003", "SIRV8001", "SIRV8002", "SIRV8003", "SIRV10001", "SIRV10002", "SIRV10003", "SIRV12001", "SIRV12002", "SIRV12003"]
    def filename = bam.name.replaceAll("${options.args}", "")

    def command = ''
    for( int i=0; i<sirvname.size(); i++ ) {
        command += "julia /opt/run.jl coverage $bam SIRVomeERCCome ${leftpos[i]} ${rightpos[i]} ${filename}.${sirvname[i]}; " 
    } 
    command
}
