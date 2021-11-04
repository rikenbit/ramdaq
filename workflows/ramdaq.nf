
// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

/*
========================================================================================
    SET UP CONFIGURATION VARIABLES
========================================================================================
*/

// Configurable variables
params.adapter = params.genome ? params.genomes[ params.genome ].adapter ?: false : false

/*
========================================================================================
    SET UP VARIABLES & VALIDATE INPUTS
========================================================================================
*/

if (params.adapter) { ch_adapter = file(params.adapter, checkIfExists: true) } else { exit 1, "Adapter file not found: ${params.adapter}" }


//Create a channel for input read files
if (params.readPaths) {
    if (params.single_end) {
        ch_reads = Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
    } else {
        ch_reads = Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
    }
} else {
    ch_reads = Channel
        .fromFilePairs(params.reads, size: params.single_end ? 1 : 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
}

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { FASTQC as FASTQC_RAW } from '../modules/local/fastqc' addParams( options: modules['fastqc'] )
include { FASTQC as FASTQC_TRIM } from '../modules/local/fastqc' addParams( options: modules['fastqc_trim'] )
include { FASTQMCF } from '../modules/local/fastqmcf' addParams( options: modules['fastqmcf'] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow RAMDAQ {
    //ch_versions = Channel.empty()

    //
    // MODULE: Read QC
    //
    FASTQC_RAW (
        ch_reads
    )
    //ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: Adapter trimming
    //
    FASTQMCF (
        ch_reads,
        params.adapter
    )
    //ch_versions = ch_versions.mix(BBMAP_BBSPLIT.out.versions.first())

    //
    // MODULE: Trimmed Read QC
    //
    FASTQC_TRIM (
        FASTQMCF.out.trimmed_reads
    )
    //ch_versions = ch_versions.mix(FASTQC.out.versions)



}

