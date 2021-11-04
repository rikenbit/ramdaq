
// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

/*
========================================================================================
    SET UP CONFIGURATION VARIABLES
========================================================================================
*/

// Configurable variables
params.adapter = params.genome ? params.genomes[ params.genome ].adapter ?: false : false
params.hisat2_idx = params.genome ? params.genomes[ params.genome ].hisat2_idx ?: false : false
params.hisat2_rrna_idx = params.genome ? params.genomes[ params.genome ].hisat2_rrna_idx ?: false : false

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

if (params.hisat2_idx) {
    if (params.hisat2_idx.endsWith('.tar.gz')) {

        file(params.hisat2_idx, checkIfExists: true)
        ch_hisat2_idx = false      

    } else {
        ch_hisat2_idx = Channel
        .from(params.hisat2_idx)
        .flatMap{file(params.hisat2_idx, checkIfExists: true)}
        .ifEmpty { exit 1, "HISAT2 index files not found: ${params.hisat2_idx}" }   
    }
}

if (params.hisat2_rrna_idx) {
    if (params.hisat2_rrna_idx.endsWith('.tar.gz')) {
        file(params.hisat2_rrna_idx, checkIfExists: true)
        ch_hisat2_rrna_idx = false
    } else {
        ch_hisat2_rrna_idx = Channel
        .from(params.hisat2_rrna_idx)
        .flatMap{file(params.hisat2_rrna_idx, checkIfExists: true)}
        .ifEmpty { exit 1, "HISAT2 rrna index files not found: ${params.hisat2_rrna_idx}" }
    }
}

/*
========================================================================================
    SET UP DIR/FILE PATH VARIABLES
========================================================================================
*/

// Tools dir
ch_tools_dir = workflow.scriptFile.parent + "/tools"


/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { UNTAR_INDEX as UNTAR_HISAT2_IDX } from '../modules/local/untar_index' addParams( options: modules['untar_index'] )
include { UNTAR_INDEX as UNTAR_HISAT2_RRNA_IDX } from '../modules/local/untar_index' addParams( options: modules['untar_index'] )
include { FASTQC as FASTQC_RAW } from '../modules/local/fastqc' addParams( options: modules['fastqc'] )
include { FASTQC as FASTQC_TRIM } from '../modules/local/fastqc' addParams( options: modules['fastqc_trim'] )
include { FASTQMCF } from '../modules/local/fastqmcf' addParams( options: modules['fastqmcf'] )
include { HISAT2 as HISAT2_ALLGENES } from '../modules/local/hisat2' addParams( options: modules['hisat2_allgenes'] )
include { HISAT2 as HISAT2_RRNA } from '../modules/local/hisat2' addParams( options: modules['hisat2_rrna'] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow RAMDAQ {
    //ch_versions = Channel.empty()

    if (!ch_hisat2_idx){
        //
        // MODULE: untar index.tar.gz [all genes]
        //
        UNTAR_HISAT2_IDX (
            params.hisat2_idx
        )
        .index_files
        .set { ch_hisat2_idx }
    }
    if (!ch_hisat2_rrna_idx){
        //
        // MODULE: untar index.tar.gz [rrna]
        //
        UNTAR_HISAT2_RRNA_IDX (
            params.hisat2_rrna_idx
        )
        .index_files
        .set { ch_hisat2_rrna_idx }
    }

    //
    // MODULE: Read QC
    //
    FASTQC_RAW (
        ch_reads
    )

    //
    // MODULE: Adapter trimming
    //
    FASTQMCF (
        ch_reads,
        params.adapter
    )
    .trimmed_reads
    .set { ch_trimmed_reads }

    //
    // MODULE: Trimmed Read QC
    //
    FASTQC_TRIM (
        ch_trimmed_reads
    )

    //
    // MODULE: Alignment with Hisat2 [all genes]
    //
    HISAT2_ALLGENES (
        ch_trimmed_reads,
        ch_hisat2_idx.collect(),
        ch_tools_dir
    )
    ch_hisat2_bam = HISAT2_ALLGENES.out.hisat2_bam_qc

    //
    // MODULE: Alignment with Hisat2 [rrna]
    //
    HISAT2_RRNA (
        ch_trimmed_reads,
        ch_hisat2_rrna_idx.collect(),
        ch_tools_dir
    )

}

