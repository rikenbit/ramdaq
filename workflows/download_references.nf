// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    SET UP CONFIGURATION VARIABLES
========================================================================================
*/

//debug
//println "  genomes '${params.genomes}' "


// Check if annot_ver exists in the config file
if (params.genomes && params.annot_ver && !params.genomes.containsKey(params.annot_ver)) {
    exit 1, "The provided annot_ver '${params.annot_ver}' is not available. Currently the available annot_ver are ${params.genomes.keySet().join(", ")}"
}

if (params.dl_references && params.outdir == './results') {
    exit 1, "Please setting the output directory name where the dl files will be saved. ex) --outdir ramdaq_dl_references"
}


// Configurable variables
params.species = params.annot_ver ? params.genomes[ params.annot_ver ].species ?: false : false

params.common = params.annot_ver ? params.genomes[ params.annot_ver ].common_files ?: false : false
params.common_md5sum = params.annot_ver ? params.genomes[ params.annot_ver ].common_files_md5sum ?: false : false
if (params.common) ch_common_dirname = params.common.replaceAll("\\.tar\\.gz", "")

params.hisat2_idx = params.annot_ver ? params.genomes[ params.annot_ver ].hisat2_idx ?: false : false
params.hisat2_idx_md5sum = params.annot_ver ? params.genomes[ params.annot_ver ].hisat2_idx_md5sum ?: false : false

params.hisat2_rrna_idx = params.annot_ver ? params.genomes[ params.annot_ver ].hisat2_rrna_idx ?: false : false
params.hisat2_rrna_idx_md5sum = params.annot_ver ? params.genomes[ params.annot_ver ].hisat2_rrna_idx_md5sum ?: false : false

params.annotation_gtf = params.annot_ver ? params.genomes[ params.annot_ver ].annotation_gtf ?: false : false
params.annotation_gtf_md5sum = params.annot_ver ? params.genomes[ params.annot_ver ].annotation_gtf_md5sum ?: false : false
if (params.annotation_gtf) annotation_gtf_dirname = params.annotation_gtf.replaceAll("\\.tar\\.gz", "")

params.hisat2_sirv_idx = params.annot_ver ? params.genomes[ params.annot_ver ].hisat2_sirv_idx ?: false : false
params.hisat2_sirv_idx_md5sum = params.annot_ver ? params.genomes[ params.annot_ver ].hisat2_sirv_idx_md5sum ?: false : false

params.rsem_allgene_idx = params.annot_ver ? params.genomes[ params.annot_ver ].rsem_allgene_idx ?: false : false
params.rsem_allgene_idx_md5sum = params.annot_ver ? params.genomes[ params.annot_ver ].rsem_allgene_idx_md5sum ?: false : false

params.rsem_sirv_idx = params.annot_ver ? params.genomes[ params.annot_ver ].rsem_sirv_idx ?: false : false
params.rsem_sirv_idx_md5sum = params.annot_ver ? params.genomes[ params.annot_ver ].rsem_sirv_idx_md5sum ?: false : false

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    ${c_green}============================================================${c_reset}
            ${c_blue}ramdaq v${workflow.manifest.version} DL_REFERENCES${c_reset} 
    ${c_green}============================================================${c_reset}
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "============================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run rikenbit/ramdaq --dl_references

    Pipeline setting:
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: docker, singularity, test, and more
      -c                              Specify the path to a specific config file
      --annot_ver [str]                  Name of human or mouse reference to download : ${params.genomes.keySet().join(", ")}
        
    Other:
      --outdir [str]                 The output directory where the dl files will be saved. ex) ramdaq_dl_references
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
      -resume                         Specify this when restarting a pipeline

    """.stripIndent()
}


/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { WGET_REFERENCES as WGET_COMMON_FILES } from '../modules/local/wget_references' 
include { UNTAR_REFERENCES as UNTAR_COMMON_FILES } from '../modules/local/untar_references' addParams( options: modules['untar_common_files'] )

include { WGET_REFERENCES as WGET_HISAT2_IDX } from '../modules/local/wget_references'
include { UNTAR_REFERENCES as UNTAR_HISAT2_IDX } from '../modules/local/untar_references' addParams( options: modules['untar_index_hisat2'] )

include { WGET_REFERENCES as WGET_HISAT2_RRNA_IDX } from '../modules/local/wget_references' 
include { UNTAR_REFERENCES as UNTAR_HISAT2_RRNA_IDX } from '../modules/local/untar_references' addParams( options: modules['untar_index_hisat2'] )

include { WGET_REFERENCES as WGET_GTF_FILES } from '../modules/local/wget_references' 
include { UNTAR_REFERENCES as UNTAR_GTF_FILES } from '../modules/local/untar_references' addParams( options: modules['untar_gtf_files'] )

include { WGET_REFERENCES as WGET_HISAT2_SIRV_IDX } from '../modules/local/wget_references' 
include { UNTAR_REFERENCES as UNTAR_HISAT2_SIRV_IDX } from '../modules/local/untar_references' addParams( options: modules['untar_index_hisat2'] )

include { WGET_REFERENCES as WGET_RSEM_ALLGENE_IDX } from '../modules/local/wget_references' 
include { UNTAR_REFERENCES as UNTAR_RSEM_ALLGENE_IDX } from '../modules/local/untar_references' addParams( options: modules['untar_index_rsem'] )

include { WGET_REFERENCES as WGET_RSEM_SIRV_IDX } from '../modules/local/wget_references' 
include { UNTAR_REFERENCES as UNTAR_RSEM_SIRV_IDX } from '../modules/local/untar_references' addParams( options: modules['untar_index_rsem'] )


/*
========================================================================================
    CREATE WORKFLOW SUMMARY
========================================================================================
*/

if (params.dl_references) {

// Header log info
log.info nfcoreHeader()

// create summary
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name'] = custom_runName ?: workflow.runName

if (params.dl_directory_url) summary['Download URL'] = params.dl_directory_url
if (params.species) summary['Species'] = params.species
summary['References ver'] = params.annot_ver

summary['Resource allocation for the entire workflow']  = "$params.entire_max_cpus cpus, $params.entire_max_memory memory"
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Save dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName

summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[32m------------------------------------------------------------\033[0m"
log.info "\033[32m------------------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'ramdaq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'ramdaq Workflow Summary'
    section_href: 'https://github.com/rikenbit/ramdaq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

}


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow DL_REFERENCES {

    def ch_empty_dirname = ""

    //
    // MODULE: untar tar.gz [common files]
    //
    WGET_COMMON_FILES (
        params.common
    )
    ch_common_md5sum_dl = WGET_COMMON_FILES.out.md5sum
    ch_common_gz = WGET_COMMON_FILES.out.dl_file

    ch_common_md5sum_dl
    .map { file(it).text.trim() }
    .ifEmpty { exit 1, "No MD5 sum was produced." }
    .filter { it == params.common_md5sum }
    .ifEmpty { exit 1, "DL did not complete successfully. : ${params.common}" }

    UNTAR_COMMON_FILES (
        ch_common_gz,
        ch_common_dirname
    )

    //
    // MODULE: untar tar.gz [hisat2 all genome index]
    //
    WGET_HISAT2_IDX (
        params.hisat2_idx
    )
    ch_hisat2_idx_md5sum_dl = WGET_HISAT2_IDX.out.md5sum
    ch_hisat2_idx_gz = WGET_HISAT2_IDX.out.dl_file

    ch_hisat2_idx_md5sum_dl
    .map { file(it).text.trim() }
    .ifEmpty { exit 1, "No MD5 sum was produced." }
    .filter { it == params.hisat2_idx_md5sum }
    .ifEmpty { exit 1, "DL did not complete successfully. : ${params.hisat2_idx}" }

    UNTAR_HISAT2_IDX (
        ch_hisat2_idx_gz,
        ch_empty_dirname
    )

    //
    // MODULE: untar tar.gz [hisat2 rrna genome index]
    //
    WGET_HISAT2_RRNA_IDX (
        params.hisat2_rrna_idx
    )
    ch_hisat2_rrna_idx_md5sum_dl = WGET_HISAT2_RRNA_IDX.out.md5sum
    ch_hisat2_rrna_idx_gz = WGET_HISAT2_RRNA_IDX.out.dl_file

    ch_hisat2_rrna_idx_md5sum_dl
    .map { file(it).text.trim() }
    .ifEmpty { exit 1, "No MD5 sum was produced." }
    .filter { it == params.hisat2_rrna_idx_md5sum }
    .ifEmpty { exit 1, "DL did not complete successfully. : ${params.hisat2_rrna_idx}" }

    UNTAR_HISAT2_RRNA_IDX (
        ch_hisat2_rrna_idx_gz,
        ch_empty_dirname
    )

    //
    // MODULE: untar tar.gz [annotation gtf files]
    //
    WGET_GTF_FILES (
        params.annotation_gtf
    )
    ch_gtf_md5sum_dl = WGET_GTF_FILES.out.md5sum
    ch_gtf_gz = WGET_GTF_FILES.out.dl_file

    ch_gtf_md5sum_dl
    .map { file(it).text.trim() }
    .ifEmpty { exit 1, "No MD5 sum was produced." }
    .filter { it == params.annotation_gtf_md5sum }
    .ifEmpty { exit 1, "DL did not complete successfully. : ${params.annotation_gtf}" }

    UNTAR_GTF_FILES (
        ch_gtf_gz,
        annotation_gtf_dirname
    )

    //
    // MODULE: untar tar.gz [hisat2 SIRV index]
    //
    WGET_HISAT2_SIRV_IDX (
        params.hisat2_sirv_idx
    )
    ch_hisat2_sirv_idx_md5sum_dl = WGET_HISAT2_SIRV_IDX.out.md5sum
    ch_hisat2_sirv_idx_gz = WGET_HISAT2_SIRV_IDX.out.dl_file

    ch_hisat2_sirv_idx_md5sum_dl
    .map { file(it).text.trim() }
    .ifEmpty { exit 1, "No MD5 sum was produced." }
    .filter { it == params.hisat2_sirv_idx_md5sum }
    .ifEmpty { exit 1, "DL did not complete successfully. : ${params.hisat2_sirv_idx}" }

    UNTAR_HISAT2_SIRV_IDX (
        ch_hisat2_sirv_idx_gz,
        ch_empty_dirname
    )

    //
    // MODULE: untar tar.gz [RSEM all genes index]
    //
    WGET_RSEM_ALLGENE_IDX (
        params.rsem_allgene_idx
    )
    ch_rsem_allgene_idx_md5sum_dl = WGET_RSEM_ALLGENE_IDX.out.md5sum
    ch_rsem_allgene_idx_gz = WGET_RSEM_ALLGENE_IDX.out.dl_file

    ch_rsem_allgene_idx_md5sum_dl
    .map { file(it).text.trim() }
    .ifEmpty { exit 1, "No MD5 sum was produced." }
    .filter { it == params.rsem_allgene_idx_md5sum }
    .ifEmpty { exit 1, "DL did not complete successfully. : ${params.rsem_allgene_idx}" }

    UNTAR_RSEM_ALLGENE_IDX (
        ch_rsem_allgene_idx_gz,
        ch_empty_dirname
    )

    //
    // MODULE: untar tar.gz [RSEM sirv index]
    //
    WGET_RSEM_SIRV_IDX (
        params.rsem_sirv_idx
    )
    ch_rsem_sirv_idx_md5sum_dl = WGET_RSEM_SIRV_IDX.out.md5sum
    ch_rsem_sirv_idx_gz = WGET_RSEM_SIRV_IDX.out.dl_file

    ch_rsem_sirv_idx_md5sum_dl
    .map { file(it).text.trim() }
    .ifEmpty { exit 1, "No MD5 sum was produced." }
    .filter { it == params.rsem_sirv_idx_md5sum }
    .ifEmpty { exit 1, "DL did not complete successfully. : ${params.rsem_sirv_idx}" }

    UNTAR_RSEM_SIRV_IDX (
        ch_rsem_sirv_idx_gz,
        ch_empty_dirname
    )

}



///////////////////////////////////////////////////////////////////////////////
/*
* workflow.onComplete
*/
///////////////////////////////////////////////////////////////////////////////

workflow.onComplete {

if (params.dl_references) {
    
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    log.info "${c_green}------------------------------------------------------------${c_reset}"
    log.info "${c_green}------------------------------------------------------------${c_reset}"

    if (workflow.success) {
        log.info "${c_purple}[ramdaq]${c_green} Pipeline completed successfully!${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[ramdaq]${c_red} Pipeline completed with errors!${c_reset}"
        log.info "${c_purple}[ramdaq]${c_red} If you get an error when using ramdaq, please refer to the foillowing:${c_reset}"
        log.info "${c_purple}[ramdaq]${c_red}    ・https://github.com/rikenbit/ramdaq/blob/master/docs/troubleshooting.md ${c_reset}"
        log.info "${c_purple}[ramdaq]${c_red}    ・https://nf-co.re/usage/troubleshooting ${c_reset}"
    }
    
    // copy .nextflow.log
    today = new Date().format("yyyy-MM-dd-HH-mm-ss")
    new File("${params.outdir}/ramdaq-${today}.log") << new File('.nextflow.log').text
    
    println "${c_purple}[ramdaq]${c_green} The log file .nextflow.log was copied to ${params.outdir}/ramdaq-${today}.log${c_reset}"
}
}

/*
========================================================================================
    THE END
========================================================================================
*/












