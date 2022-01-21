
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

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Configurable variables
params.adapter = params.genome ? params.genomes[ params.genome ].adapter ?: false : false
params.hisat2_idx = params.genome ? params.genomes[ params.genome ].hisat2_idx ?: false : false
params.hisat2_rrna_idx = params.genome ? params.genomes[ params.genome ].hisat2_rrna_idx ?: false : false
params.chrsize = params.genome ? params.genomes[ params.genome ].chrsize ?: false : false
params.bed = params.genome ? params.genomes[ params.genome ].bed ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.mt_gtf = params.genome ? params.genomes[ params.genome ].mt_gtf ?: false : false
params.histone_gtf = params.genome ? params.genomes[ params.genome ].histone_gtf ?: false : false
params.rsem_allgene_idx = params.genome ? params.genomes[ params.genome ].rsem_allgene_idx ?: false : false

// option: SIRV quantification
params.hisat2_sirv_idx = params.genome ? params.genomes[ params.genome ].hisat2_sirv_idx ?: false : false
params.rsem_sirv_idx = params.genome ? params.genomes[ params.genome ].rsem_sirv_idx ?: false : false

/*
========================================================================================
    SET UP VARIABLES & VALIDATE INPUTS
========================================================================================
*/

if (params.stranded && params.stranded != 'unstranded' && params.stranded != 'fr-firststrand' && params.stranded != 'fr-secondstrand') {
    exit 1, "Invalid stranded option: ${params.stranded}. Valid options: 'unstranded' or 'fr-firststrand' or 'fr-secondstrand'!"
}

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

if (params.chrsize) { ch_chrsize= file(params.chrsize, checkIfExists: true) } else { exit 1, "Chromosome sizes file not found: ${params.chrsize}" }
if (params.bed) { ch_bed= file(params.bed, checkIfExists: true) } else { exit 1, "BED file not found: ${params.bed}" }
if (params.gtf) { ch_gtf= file(params.gtf, checkIfExists: true) } else { exit 1, "GTF annotation file not found: ${params.gtf}" }
if (params.mt_gtf) { ch_mt_gtf= file(params.mt_gtf, checkIfExists: true) } else { exit 1, "Mitocondria GTF annotation file not found: ${params.mt_gtf}" }
if (params.histone_gtf) { ch_histone_gtf= file(params.histone_gtf, checkIfExists: true) } else { exit 1, "Histone GTF annotation file not found: ${params.histone_gtf}" }

if (params.rsem_allgene_idx) {
    if (params.rsem_allgene_idx.endsWith('.tar.gz')) {

        file(params.rsem_allgene_idx, checkIfExists: true)
        ch_rsem_allgene_idx = false      

    } else {
        ch_rsem_allgene_idx = Channel
        .from(params.rsem_allgene_idx)
        .flatMap{file(params.rsem_allgene_idx, checkIfExists: true)}
        .ifEmpty { exit 1, "RSEM All genes files not found: ${params.rsem_allgene_idx}" }   
    }
}

// option: ERCC / SIRV quantification
if (params.spike_in_ercc && params.spike_in_ercc.toString() == 'true') { ch_spike_in_ercc = params.spike_in_ercc_default_amount } else { ch_spike_in_ercc = params.spike_in_ercc }
if (params.spike_in_sirv && params.spike_in_sirv.toString() == 'true') { exit 1, "--spike_in_sirv option requires a dilution rate value (e.g. --spike_in_sirv '4e-6')" }
if (params.spike_in_sirv) { ch_spike_in_ercc = params.spike_in_sirv }

if (params.spike_in_sirv) {

    if (params.hisat2_sirv_idx) {
        if (params.hisat2_sirv_idx.endsWith('.tar.gz')) {
            file(params.hisat2_sirv_idx, checkIfExists: true)
            ch_hisat2_sirv_idx = false      

        } else {
            ch_hisat2_sirv_idx = Channel
            .from(params.hisat2_sirv_idx)
            .flatMap{file(params.hisat2_sirv_idx, checkIfExists: true)}
            .ifEmpty { exit 1, "HISAT2 SIRVome index files not found: ${params.hisat2_sirv_idx}" }   
        }
    }

    if (params.rsem_sirv_idx) {
        if (params.rsem_sirv_idx.endsWith('.tar.gz')) {

            file(params.rsem_sirv_idx, checkIfExists: true)
            ch_rsem_sirv_idx = false      

        } else {
            ch_rsem_sirv_idx = Channel
            .from(params.rsem_sirv_idx)
            .flatMap{file(params.rsem_sirv_idx, checkIfExists: true)}
            .ifEmpty { exit 1, "RSEM SIRVome index files not found: ${params.rsem_sirv_idx}" }   
        }
    }

} else {
    ch_hisat2_sirv_idx = false
    ch_rsem_sirv_idx = false
}


/*
========================================================================================
    SET UP HEADERS FOR REPORT
========================================================================================
*/

ch_biotypes_header = file("$projectDir/assets/biotypes_header.txt", checkIfExists: true)
ch_mdsplot_header = file("$projectDir/assets/mdsplot_header.txt", checkIfExists: true)
ch_heatmap_header = file("$projectDir/assets/heatmap_header.txt", checkIfExists: true)
ch_ercc_data = params.spike_in_sirv ? file("$projectDir/assets/ercc_in_sirv_dataset.txt", checkIfExists: true) : file("$projectDir/assets/ercc_dataset.txt", checkIfExists: true)
ch_ercc_corr_header = file("$projectDir/assets/ercc_correlation_header.txt", checkIfExists: true)
ch_ercc_corr_header_gstat = file("$projectDir/assets/gstat_ercc_correlation_header.txt", checkIfExists: true)
ch_assignedgenome_header = file("$projectDir/assets/barplot_assignedgenome_rate_header.txt", checkIfExists: true)
ch_assignedgenome_header_gstat = file("$projectDir/assets/gstat_assignedgenome_rate_header.txt", checkIfExists: true)

ch_fcounts_allgene_header = file("$projectDir/assets/barplot_fcounts_allgene_header.txt", checkIfExists: true)
ch_fcounts_allgene_header_gstat = file("$projectDir/assets/gstat_fcounts_allgene_header.txt", checkIfExists: true)
ch_fcounts_mt_header = file("$projectDir/assets/barplot_fcounts_mt_header.txt", checkIfExists: true)
ch_fcounts_mt_header_gstat = file("$projectDir/assets/gstat_fcounts_mt_header.txt", checkIfExists: true)
ch_fcounts_histone_header = file("$projectDir/assets/barplot_fcounts_histone_header.txt", checkIfExists: true)
ch_fcounts_histone_header_gstat = file("$projectDir/assets/gstat_fcounts_histone_header.txt", checkIfExists: true)

ch_num_of_detgene_header = file("$projectDir/assets/barplot_num_of_detgene_header.txt", checkIfExists: true)
ch_num_of_detgene_header_gstat = file("$projectDir/assets/gstat_num_of_detgene_header.txt", checkIfExists: true)
ch_pcaplot_header = file("$projectDir/assets/pcaplot_header.txt", checkIfExists: true)
ch_tsneplot_header = file("$projectDir/assets/tsneplot_header.txt", checkIfExists: true)
ch_umapplot_header = file("$projectDir/assets/umapplot_header.txt", checkIfExists: true)

ch_num_of_gene_rsem_header = file("$projectDir/assets/barplot_num_of_gene_rsem_header.txt", checkIfExists: true)
ch_num_of_gene_rsem_header_gstat = file("$projectDir/assets/gstat_num_of_gene_rsem_header.txt", checkIfExists: true)
ch_num_of_ts_rsem_header = file("$projectDir/assets/barplot_num_of_ts_rsem_header.txt", checkIfExists: true)
ch_num_of_ts_rsem_header_gstat = file("$projectDir/assets/gstat_num_of_ts_rsem_header.txt", checkIfExists: true)

ch_entropy_of_sirv_header = Channel.fromPath("$projectDir/assets/barplot_entropy_of_sirv_header.txt", checkIfExists: true)
ch_entropy_of_sirv_header_gstat = Channel.fromPath("$projectDir/assets/gstat_entropy_of_sirv_header.txt", checkIfExists: true)

/*
========================================================================================
    SET UP DIR/FILE PATH VARIABLES
========================================================================================
*/

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : Channel.empty()

// Tools dir
ch_tools_dir = workflow.scriptFile.parent + "/tools"

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

// Function that checks the number of mapped reads from flagstat output
// and returns true if > params.min_mapped_reads and otherwise false
def check_mappedread(bam,flagstat_path,min_mapped_reads=10) {

    def filename = bam.getName()
    def flagstat_file = new File(flagstat_path.toString())
    def mapped = 0
    flagstat_file.eachLine { line ->
        if (line.contains(' mapped (')) {
            mapped = line.tokenize().first().toInteger()
        }
    }
    if (mapped < min_mapped_reads.toInteger()) {
        log.warn "[ramdaq] ${filename} Failed mapped reads threshold: ${mapped} < ${params.min_mapped_reads}. Ignoring downstream analysis."
        return false
    } else {
        return true
    }
}

def check_mappedread_sirv(bam,flagstat_path,min_mapped_reads=10) {
    
    def filename = bam.getName()

    //trim R1 or R2 bams
    if (filename.indexOf(".R1.") > 0 || filename.indexOf(".R2.") > 0){
        return false
    }
    
    def flagstat_file = new File(flagstat_path.toString())
    def mapped = 0
    flagstat_file.eachLine { line ->
        if (line.contains(' mapped (')) {
            mapped = line.tokenize().first().toInteger()
        }
    }
    if (mapped < min_mapped_reads.toInteger()) {
        log.warn "[ramdaq] ${filename} Failed mapped reads threshold: ${mapped} < ${params.min_mapped_reads}. Ignoring downstream SIRV quantification."
        return false
    } else {
        return true
    }
}

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
            ${c_purple}ramdaq v${workflow.manifest.version}${c_reset}
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

    nextflow run rikenbit/ramdaq --reads '*_R{1,2}.fastq.gz' -profile docker

    Pipeline setting:
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: docker, singularity, test, and more
      -c                              Specify the path to a specific config file
      --reads [file]                  Path to input data (must be surrounded with quotes)
      --single_end                    Specifies that the input is single-end reads
      --stranded [str]                unstranded : default
                                      fr-firststrand : First read corresponds to the reverse complemented counterpart of a transcript
                                      fr-secondstrand : First read corresponds to a transcript
      --genome [str]                  Name of human or mouse latest reference: ${params.genomes.keySet().join(", ")}
      --saveReference                 Save the generated reference files to the results directory
      --local_annot_dir [str]         Base path for local annotation files
      --entire_max_cpus [N]            Maximum number of CPUs to use for each step of the pipeline. Should be in form e.g. --entire_max_cpus 16. Default: '${params.entire_max_cpus}'
      --entire_max_memory [str]        Memory limit for each step of pipeline. Should be in form e.g. --entire_max_memory '16.GB'. Default: '${params.entire_max_memory}'  
        
    Other:
      --outdir [str]                  The output directory where the results will be saved
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
      -resume                         Specify this when restarting a pipeline
      --max_memory [str]              Memory limit for each step of pipeline. Should be in form e.g. --max_memory '8.GB'. Default: '${params.max_memory}'
      --max_time [str]                Time limit for each step of the pipeline. Should be in form e.g. --max_time '2.h'. Default: '${params.max_time}'
      --max_cpus [str]                Maximum number of CPUs to use for each step of the pipeline. Should be in form e.g. --max_cpus 1. Default: '${params.max_cpus}'
      --monochrome_logs               Set to disable colourful command line output and live life in monochrome

    Fastqmcf:
      --maxReadLength [N]             Maximum remaining sequence length (Default: 75)
      --minReadLength [N]             Minimum remaining sequence length (Default: 36)
      --skew [N]                      Skew percentage-less-than causing cycle removal (Default: 4)
      --quality [N]                   Quality threshold causing base removal (Default: 30)
    
    Hisat2:
      --softclipping                  HISAT2 allow soft-clip reads near their 5' and 3' ends (Default: disallow)
      --hs_threads_num [N]            HISAT2 to launch a specified number of parallel search threads (Default: 1)
    
    RSEM:
      --rsem_threads_num [N]          Number of threads to use (Default: 1)
    
    FeatureCounts:
      --extra_attributes              Define which extra parameters should also be included in featureCounts (Default: 'gene_name')
      --group_features                Define the attribute type used to group features (Default: 'gene_id')
      --count_type                    Define the type used to assign reads (Default: 'exon')
      --allow_multimap                Multi-mapping reads/fragments will be counted (Default: true)
      --allow_overlap                 Reads will be allowed to be assigned to more than one matched meta-feature (Default: true)
      --count_fractionally            Assign fractional counts to features  (Default: true / This option must be used together with ‘--allow_multimap’ or ‘--allow_overlap’ or both)
      --fc_threads_num [N]            Number of the threads (Default: 1)
      --group_features_type           Define the type attribute used to group features based on the group attribute (default: 'gene_type')
    
    For ERCC RNA Spike-In Controls:
      --spike_in_ercc [str]           Dilution rate of the ERCC Spike-In Control Mix 1. Use when the samples contain the ERCC Spike-In Control Mix 1. The value is used to calculate copy number of ERCC. If the value is not specified, '2e-7' is used as dilution rate. (default: false)
      --spike_in_sirv [str]           Dilution rate of the SIRV-Set 4. Use when the samples contain the SIRV-Set 4. The value is used to calculate copy number of ERCC in the SIRV-Set 4. (default: false)
    
    MultiQC report:
      --sampleLevel                   Used to turn off the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples

    """.stripIndent()
}


/*
========================================================================================
    CREATE WORKFLOW SUMMARY
========================================================================================
*/

// Header log info
log.info nfcoreHeader()

// create summary
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name'] = custom_runName ?: workflow.runName
if (!params.readPaths) summary['Reads'] = params.reads
summary['Data Type'] = params.single_end ? 'Single-End' : 'Paired-End'
if (params.stranded)  {
    if (params.stranded == 'unstranded') summary['Strandness'] = 'Unstranded'
    if (params.stranded == 'fr-firststrand') summary['Strandness'] = 'Forward stranded'
    if (params.stranded == 'fr-secondstrand') summary['Strandness'] = 'Reverse stranded'
} else {
    summary['Strandness'] = 'Unstranded'
}
summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
if (params.hisat2_idx) summary['HISAT2 Index'] = params.hisat2_idx
if (params.hisat2_rrna_idx) summary['HISAT2 rRNA Index'] = params.hisat2_rrna_idx
if (params.rsem_allgene_idx) summary['RSEM-Bowtie2 All genes Index'] = params.rsem_allgene_idx
if (params.chrsize)  summary['Chromosome sizes'] = params.chrsize
if (params.bed) summary['BED Annotation'] = params.bed
if (params.gtf) summary['GTF Annotation'] = params.gtf
if (params.mt_gtf) summary['Mitocondria GTF Annotation'] = params.mt_gtf
if (params.histone_gtf) summary['Histone GTF Annotation'] = params.histone_gtf
// featureCounts options
if (params.allow_multimap) summary['Multimap Reads'] = params.allow_multimap ? 'Allow' : 'Disallow'
if (params.allow_overlap) summary['Overlap Reads'] = params.allow_overlap ? 'Allow' : 'Disallow'
if (params.count_fractionally) summary['Fractional counting'] = params.count_fractionally ? 'Enabled' : 'Disabled'
if (params.group_features_type) summary['Biotype GTF field'] = params.group_features_type

if (params.min_mapped_reads) summary['Min Mapped Reads'] = params.min_mapped_reads
summary['ERCC quantification mode']   = params.spike_in_ercc || params.spike_in_sirv ? 'On' : 'Off'
summary['SIRV quantification mode']   = params.spike_in_sirv ? 'On' : 'Off'
if (params.hisat2_sirv_idx) summary['HISAT2 SIRVome Index'] = params.hisat2_sirv_idx
if (params.rsem_sirv_idx) summary['RSEM-Bowtie2 SIRVome Index'] = params.rsem_sirv_idx

summary['Resource allocation for the entire workflow']  = "$params.entire_max_cpus cpus, $params.entire_max_memory memory"
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
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
log.info "\033[32m----------------------------------------------------\033[0m"
log.info "\033[32m----------------------------------------------------\033[0m"

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

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: modules['get_software_versions'] )

include { UNTAR_INDEX as UNTAR_HISAT2_IDX } from '../modules/local/untar_index' addParams( options: modules['untar_index_hisat2'] )
include { UNTAR_INDEX as UNTAR_HISAT2_RRNA_IDX } from '../modules/local/untar_index' addParams( options: modules['untar_index_hisat2'] )
include { UNTAR_INDEX as UNTAR_RSEM_ALLGENE_IDX } from '../modules/local/untar_index' addParams( options: modules['untar_index_rsem'] )

include { UNTAR_INDEX as UNTAR_HISAT2_SIRV_IDX } from '../modules/local/untar_index' addParams( options: modules['untar_index_hisat2'] )
include { UNTAR_INDEX as UNTAR_RSEM_SIRV_IDX } from '../modules/local/untar_index' addParams( options: modules['untar_index_rsem'] )

include { FASTQC as FASTQC_RAW } from '../modules/local/fastqc' addParams( options: modules['fastqc'] )
include { FASTQC as FASTQC_TRIM } from '../modules/local/fastqc' addParams( options: modules['fastqc_trim'] )
include { FASTQMCF } from '../modules/local/fastqmcf' addParams( options: modules['fastqmcf'] )

include { HISAT2 as HISAT2_GENOME } from '../modules/local/hisat2' addParams( options: modules['hisat2_genome'] )
include { HISAT2 as HISAT2_RRNA } from '../modules/local/hisat2' addParams( options: modules['hisat2_rrna'] )
include { MERGE_SUMMARYFILE as MERGE_SUMMARYFILE_HISAT2 } from '../modules/local/merge_summaryfile' addParams( options: modules['merge_summaryfile_hisat2'] )

include { BAM2WIG as BAM2WIG_ALLGENES } from '../modules/local/bam2wig' addParams( options: modules['bam2wig'] )
include { ADJUST_BED_NONCODING } from '../modules/local/adjust_bed_noncoding' addParams( options: modules['adjust_bed_noncoding'] )
include { RSEQC } from '../modules/local/rseqc' addParams( options: modules['rseqc'] )
include { MERGE_SUMMARYFILE as MERGE_SUMMARYFILE_READDIST } from '../modules/local/merge_summaryfile' addParams( options: modules['merge_summaryfile_readdist'] )
include { READCOVERAGE } from '../modules/local/readcoverage' addParams( options: modules['readcoverage'] )

include { FEATURECOUNTS as FEATURECOUNTS_ALL_GTF } from '../modules/local/featurecounts' addParams( options: modules['featurecounts_all_gtf'] )
include { FEATURECOUNTS as FEATURECOUNTS_MT_GTF } from '../modules/local/featurecounts' addParams( options: modules['featurecounts_mt_gtf'] )
include { FEATURECOUNTS as FEATURECOUNTS_HISTONE_GTF } from '../modules/local/featurecounts' addParams( options: modules['featurecounts_histone_gtf'] )
include { MERGE_FEATURECOUNTS as MERGE_FEATURECOUNTS_ALLGENE } from '../modules/local/merge_featurecounts' addParams( options: modules['merge_featurecounts_allgene'] )
include { MERGE_FEATURECOUNTS as MERGE_FEATURECOUNTS_MT } from '../modules/local/merge_featurecounts' addParams( options: modules['merge_featurecounts_mt'] )
include { MERGE_FEATURECOUNTS as MERGE_FEATURECOUNTS_HISTONE } from '../modules/local/merge_featurecounts' addParams( options: modules['merge_featurecounts_histone'] )

include { RSEM_BOWTIE2 as RSEM_BOWTIE2_ALLGENES } from '../modules/local/rsem_bowtie2' addParams( options: modules['rsem_bowtie2_allgenes'] )
include { MERGE_RSEM as MERGE_RSEM_GENES } from '../modules/local/merge_rsem' addParams( options: modules['merge_rsem_genes'] )
include { MERGE_RSEM as MERGE_RSEM_ISOFORMS } from '../modules/local/merge_rsem' addParams( options: modules['merge_rsem_isoforms'] )

include { CALC_SAMPLE_CORRELATION } from '../modules/local/calc_sample_correlation' addParams( options: modules['calc_sample_correlation'] )
include { CALC_ERCC_CORRELATION } from '../modules/local/calc_ercc_correlation' addParams( options: modules['calc_ercc_correlation'] )
include { CALC_ASSIGNEDGENOME_RATE } from '../modules/local/calc_assignedgenome_rate' addParams( options: modules['calc_assignedgenome_rate'] )
include { CALC_FEATURECOUNTS_MAPRATE as CALC_FEATURECOUNTS_MAPRATE_ALLGENE } from '../modules/local/calc_featurecounts_maprate' addParams( options: modules['calc_featurecounts_maprate_allgene'] )
include { CALC_FEATURECOUNTS_MAPRATE as CALC_FEATURECOUNTS_MAPRATE_MT} from '../modules/local/calc_featurecounts_maprate' addParams( options: modules['calc_featurecounts_maprate_mt'] )
include { CALC_FEATURECOUNTS_MAPRATE as CALC_FEATURECOUNTS_MAPRATE_HISTONE} from '../modules/local/calc_featurecounts_maprate' addParams( options: modules['calc_featurecounts_maprate_histone'] )
include { CALC_DETECTEDGENES_DR as CALC_TPMCOUNTS_FEATURECOUNTS } from '../modules/local/calc_detectedgenes_dr' addParams( options: modules['calc_tpmcounts_featurecounts'] )
include { CALC_DETECTEDGENES_DR as CALC_TPMCOUNTS_RSEM_GENE } from '../modules/local/calc_detectedgenes_dr' addParams( options: modules['calc_tpmcounts_rsem_gene'] )
include { CALC_DETECTEDGENES_DR as CALC_TPMCOUNTS_RSEM_TS } from '../modules/local/calc_detectedgenes_dr' addParams( options: modules['calc_tpmcounts_rsem_ts'] )

// option: ERCC / SIRV quantification
include { HISAT2 as HISAT2_SIRV } from '../modules/local/hisat2' addParams( options: modules['hisat2_sirv'] )
include { READCOVERAGE_SIRV } from '../modules/local/readcoverage_sirv' addParams( options: modules['readcoverage_sirv'] )
include { MERGE_READCOVERAGE_SIRV } from '../modules/local/merge_readcoverage_sirv' addParams( options: modules['merge_readcoverage_sirv'] )
include { RSEM_BOWTIE2 as RSEM_BOWTIE2_SIRV } from '../modules/local/rsem_bowtie2' addParams( options: modules['rsem_bowtie2_sirv'] )
include { MERGE_RSEM as MERGE_RSEM_ISOFORMS_SIRV } from '../modules/local/merge_rsem' addParams( options: modules['merge_rsem_isoforms_sirv'] )
include { CALC_ENTROPY_SIRV } from '../modules/local/calc_entropy_sirv' addParams( options: modules['calc_entropy_sirv'] )

include { MULTIQC } from '../modules/local/multiqc' addParams( options: modules['multiqc'] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow RAMDAQ {
    //ch_versions = Channel.empty()

    //
    // MODULE: Software version output
    //
    GET_SOFTWARE_VERSIONS ()
    .software_versions_yaml
    .set { ch_software_versions_yaml }

    if (!ch_hisat2_idx){
        //
        // MODULE: untar index.tar.gz [hisat2 all genome]
        //
        UNTAR_HISAT2_IDX (
            params.hisat2_idx
        )
        .index_files
        .set { ch_hisat2_idx }
    }
    if (!ch_hisat2_rrna_idx){
        //
        // MODULE: untar index.tar.gz [hisat2 rrna genome]
        //
        UNTAR_HISAT2_RRNA_IDX (
            params.hisat2_rrna_idx
        )
        .index_files
        .set { ch_hisat2_rrna_idx }
    }
    if (!ch_rsem_allgene_idx){
        //
        // MODULE: untar index.tar.gz [rsem all genes]
        //
        UNTAR_RSEM_ALLGENE_IDX (
            params.rsem_allgene_idx
        )
        .index_files
        .set { ch_rsem_allgene_idx }
    }

    if (params.spike_in_sirv && !ch_hisat2_sirv_idx){
        //
        // MODULE: untar index.tar.gz [hisat2 sirv genome]
        //
        UNTAR_HISAT2_SIRV_IDX (
            params.hisat2_sirv_idx
        )
        .index_files
        .set { ch_hisat2_sirv_idx }
    }

    if (params.spike_in_sirv && !ch_rsem_sirv_idx){
        //
        // MODULE: untar index.tar.gz [rsem sirv genes]
        //
        UNTAR_RSEM_SIRV_IDX (
            params.rsem_sirv_idx
        )
        .index_files
        .set { ch_rsem_sirv_idx }
    }

    //
    // MODULE: Read QC
    //
    FASTQC_RAW (
        ch_reads
    ).fastqc_results
    .set { ch_fastqc_results_raw }

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
    ).fastqc_results
    .set { ch_fastqc_results_trimmed }

    //
    // MODULE: Alignment with Hisat2 [all genes]
    //
    HISAT2_GENOME (
        ch_trimmed_reads,
        ch_hisat2_idx.collect(),
        ch_tools_dir
    )
    ch_hisat2_bam_qc = HISAT2_GENOME.out.hisat2_bam_qc
    ch_hisat2_bam_count = HISAT2_GENOME.out.hisat2_bam_count
    //ch_hisat2_bam_samplenum = HISAT2_GENOME.out.hisat2_bam_samplenum
    ch_hisat2_summary = HISAT2_GENOME.out.hisat2_summary

    //
    // FUNCTION: Filtering low mapped reads bams [for bamQC, RSEM]
    //
    ch_hisat2_bam_qc
    .transpose()
    .filter { name, bam, bai, flagstat -> check_mappedread(bam,flagstat,params.min_mapped_reads) }
    .map { it[0..2] }
    .set{
        ch_hisat2_bam_qc_filtered
    }

    //
    // FUNCTION: Filtering low mapped reads bams [for featureCounts]
    //
    ch_hisat2_bam_count
    .filter { name, bam, bai, flagstat -> check_mappedread(bam,flagstat,params.min_mapped_reads) }
    .map { it[0..2] }
    .set { 
        ch_hisat2_bam_featurecount
    }

    //debug
    //ch_hisat2_bam_qc_filtered.subscribe {  println "hisat2_bam1: $it"  }
    //ch_hisat2_bam_featurecount.subscribe {  println "hisat2_bam2: $it"  }

    //
    // FUNCTION: Count the number of valid samples [for correlation plot]
    //
    ch_hisat2_bam_count
    .filter { name, bam, bai, flagstat -> check_mappedread(bam,flagstat,params.min_mapped_reads) }
    .map { it[1] }
    .set { ch_num_of_bam }
    ch_num_of_bam_list = ch_num_of_bam.toList()

    //
    // MODULE: Merge summaryfiles [Hisat2 totalseq]
    //
    MERGE_SUMMARYFILE_HISAT2 (
        ch_hisat2_summary.collect()
    )
    .merged_summary
    .set{
        ch_hisat2_merged_totalseq
    }

    //
    // MODULE: Alignment with Hisat2 [rrna]
    //
    HISAT2_RRNA (
        ch_trimmed_reads,
        ch_hisat2_rrna_idx.collect(),
        ch_tools_dir
    ).hisat2_summary
    .set { ch_hisat2_summary_rrna }

    //
    // MODULE: Bam to BigWig
    //
    BAM2WIG_ALLGENES (
        ch_hisat2_bam_qc_filtered,
        ch_chrsize
    )

    //
    // MODULE: Preparation of RseQC
    //
    ADJUST_BED_NONCODING (
        ch_bed
    )
    .bed_adjusted
    .set { ch_bed_adjusted }

    //
    // MODULE: RseQC (Bam QC)
    //
    RSEQC (
        ch_hisat2_bam_qc_filtered,
        ch_bed_adjusted
    )
    ch_readdist_totalread = RSEQC.out.readdist_totalread
    ch_rseqc_results = RSEQC.out.rseqc_results

    //
    // MODULE: Merge summaryfiles [ReadDist totalread]
    //
    MERGE_SUMMARYFILE_READDIST (
        ch_readdist_totalread.collect()
    )
    .merged_summary
    .set{
        ch_readdist_merged_totalread
    }

    //
    // MODULE: readcoverage.jl
    //
    READCOVERAGE (
        ch_hisat2_bam_qc_filtered,
        ch_bed
    ).readcov_results
    .set { ch_readcov_results }

    ch_rseqc_results_merge = ch_rseqc_results.concat(ch_readcov_results)

    //
    // MODULE: featureCounts (All-genes GTF)
    //
    FEATURECOUNTS_ALL_GTF (
        ch_hisat2_bam_featurecount,
        ch_gtf,
        ch_biotypes_header
    )
    ch_counts_to_merge_all = FEATURECOUNTS_ALL_GTF.out.counts_to_merge
    ch_counts_to_merge_all_list = ch_counts_to_merge_all.toList()
    ch_counts_summary_all = FEATURECOUNTS_ALL_GTF.out.counts_summary
    ch_counts_to_plot_corr = FEATURECOUNTS_ALL_GTF.out.counts_to_plot_corr
    ch_counts_biotype = FEATURECOUNTS_ALL_GTF.out.counts_biotype

    //
    // MODULE: Merge featureCounts output (All-genes)
    //
    MERGE_FEATURECOUNTS_ALLGENE (
        ch_counts_to_merge_all.collect()
    )
    ch_featurecounts_merged_allgene = MERGE_FEATURECOUNTS_ALLGENE.out.merged_counts
    ch_featurecounts_tpm_merged = MERGE_FEATURECOUNTS_ALLGENE.out.counts_tpm_merged
    ch_ercc_tpm_merged = MERGE_FEATURECOUNTS_ALLGENE.out.ercc_tpm_merged
    ch_ercc_tpm_merged_list = ch_ercc_tpm_merged.toList()

    //
    // MODULE: featureCounts (Mitocondria GTF)
    //
    FEATURECOUNTS_MT_GTF (
        ch_hisat2_bam_featurecount,
        ch_mt_gtf,
        ch_biotypes_header
    )
    ch_counts_to_merge_mt = FEATURECOUNTS_MT_GTF.out.counts_to_merge
    ch_counts_summary_mt = FEATURECOUNTS_MT_GTF.out.counts_summary

    //
    // MODULE: Merge featureCounts output (Mitocondria-genes)
    //
    MERGE_FEATURECOUNTS_MT (
        ch_counts_to_merge_mt.collect()
    )
    .merged_counts
    .set{
        ch_featurecounts_merged_mt
    }

    //
    // MODULE: featureCounts (Histone GTF)
    //
    FEATURECOUNTS_HISTONE_GTF (
        ch_hisat2_bam_featurecount,
        ch_histone_gtf,
        ch_biotypes_header
    )
    ch_counts_to_merge_histone = FEATURECOUNTS_HISTONE_GTF.out.counts_to_merge
    ch_counts_summary_histone = FEATURECOUNTS_HISTONE_GTF.out.counts_summary

    //
    // MODULE: Merge featureCounts output (Histone-genes)
    //
    MERGE_FEATURECOUNTS_HISTONE (
        ch_counts_to_merge_histone.collect()
    )
    .merged_counts
    .set{
        ch_featurecounts_merged_histone
    }

    //
    // MODULE: Quantification with RSEM [all genes]
    //
    RSEM_BOWTIE2_ALLGENES (
        ch_hisat2_bam_qc_filtered,
        ch_trimmed_reads,
        ch_rsem_allgene_idx.collect()
    )
    ch_rsem_isoforms_to_merge = RSEM_BOWTIE2_ALLGENES.out.rsem_isoforms_to_merge
    ch_rsem_genes_to_merge = RSEM_BOWTIE2_ALLGENES.out.rsem_genes_to_merge
    ch_rsem_results_stat = RSEM_BOWTIE2_ALLGENES.out.rsem_results_stat

    //
    // MODULE: Merge RSEM output (genes)
    //
    MERGE_RSEM_GENES (
        ch_rsem_genes_to_merge.collect()
    )
    .merged_counts
    .set{
        ch_rsem_merged_gene
    }

    //
    // MODULE: Merge RSEM output (isoforms)
    //
    MERGE_RSEM_ISOFORMS (
        ch_rsem_isoforms_to_merge.collect()
    )
    .merged_counts
    .set{
        ch_rsem_merged_ts
    }
    
    //debug
    //ch_num_of_bam_list.subscribe {  println "ch_num_of_bam: $it"  }
    //ch_num_of_bam_list.size().subscribe {  println "ch_bam_size_chk: $it"  }
    if (!params.sampleLevel) {
        //
        // MODULE: Calc sample corr
        //
        CALC_SAMPLE_CORRELATION (
            ch_counts_to_plot_corr.collect(),
            ch_num_of_bam_list.size(),
            ch_mdsplot_header,
            ch_heatmap_header
        ).sample_correlation
        .set { ch_sample_correlation }
    }

    def ch_ercc_correlation_barplot  =  Channel.empty()
    def ch_ercc_correlation_gstat  =  Channel.empty()
    if (params.spike_in_ercc || params.spike_in_sirv) {
        //
        // MODULE: Calc ERCC mol vs exp corr
        //
        CALC_ERCC_CORRELATION (
            ch_spike_in_ercc,
            ch_ercc_tpm_merged_list.size(),
            ch_ercc_tpm_merged,
            ch_ercc_data,
            ch_ercc_corr_header,
            ch_ercc_corr_header_gstat
        )
        ch_ercc_correlation_barplot = CALC_ERCC_CORRELATION.out.ercc_correlation_barplot
        ch_ercc_correlation_gstat = CALC_ERCC_CORRELATION.out.ercc_correlation_gstat
    }

    //
    // MODULE: Calc assignedgenome rate
    //
    CALC_ASSIGNEDGENOME_RATE (
        ch_hisat2_merged_totalseq,
        ch_readdist_merged_totalread,
        ch_assignedgenome_header,
        ch_assignedgenome_header_gstat
    )
    ch_assignedgenome_rate_barplot = CALC_ASSIGNEDGENOME_RATE.out.assignedgenome_rate_barplot
    ch_assignedgenome_rate_gstat = CALC_ASSIGNEDGENOME_RATE.out.assignedgenome_rate_gstat

    //
    // MODULE: Calc featureCounts mapped rate (all genes)
    //
    CALC_FEATURECOUNTS_MAPRATE_ALLGENE (
        ch_hisat2_merged_totalseq,
        ch_featurecounts_merged_allgene,
        ch_fcounts_allgene_header,
        ch_fcounts_allgene_header_gstat
    )
    ch_fcounts_maprate_barplot_allgene = CALC_FEATURECOUNTS_MAPRATE_ALLGENE.out.fcounts_maprate_barplot
    ch_fcounts_maprate_gstat_allgene = CALC_FEATURECOUNTS_MAPRATE_ALLGENE.out.fcounts_maprate_gstat

    //
    // MODULE: Calc featureCounts mapped rate (mt genes)
    //
    CALC_FEATURECOUNTS_MAPRATE_MT (
        ch_hisat2_merged_totalseq,
        ch_featurecounts_merged_mt,
        ch_fcounts_mt_header,
        ch_fcounts_mt_header_gstat
    )
    ch_fcounts_maprate_barplot_mt = CALC_FEATURECOUNTS_MAPRATE_MT.out.fcounts_maprate_barplot
    ch_fcounts_maprate_gstat_mt = CALC_FEATURECOUNTS_MAPRATE_MT.out.fcounts_maprate_gstat

    //
    // MODULE: Calc featureCounts mapped rate (hisotne genes)
    //
    CALC_FEATURECOUNTS_MAPRATE_HISTONE (
        ch_hisat2_merged_totalseq,
        ch_featurecounts_merged_histone,
        ch_fcounts_histone_header,
        ch_fcounts_histone_header_gstat
    )
    ch_fcounts_maprate_barplot_histone = CALC_FEATURECOUNTS_MAPRATE_HISTONE.out.fcounts_maprate_barplot
    ch_fcounts_maprate_gstat_histone = CALC_FEATURECOUNTS_MAPRATE_HISTONE.out.fcounts_maprate_gstat

    //
    // MODULE: Calc detected genes and dimentional reduction from featureCounts TPM
    //
    CALC_TPMCOUNTS_FEATURECOUNTS (
        ch_featurecounts_tpm_merged,
        ch_num_of_detgene_header,
        ch_num_of_detgene_header_gstat,
        ch_pcaplot_header,
        ch_tsneplot_header,
        ch_umapplot_header
    )
    ch_drplot = CALC_TPMCOUNTS_FEATURECOUNTS.out.drplot
    ch_detectedgene_barplot = CALC_TPMCOUNTS_FEATURECOUNTS.out.detectedgene_barplot
    ch_detectedgene_gstat = CALC_TPMCOUNTS_FEATURECOUNTS.out.detectedgene_gstat

    //
    // MODULE: Calc detected genes from RSEM TPM
    //
    CALC_TPMCOUNTS_RSEM_GENE (
        ch_rsem_merged_gene,
        ch_num_of_gene_rsem_header,
        ch_num_of_gene_rsem_header_gstat,
        [],
        [],
        []
    )
    ch_detectedgene_barplot_rsem_gene = CALC_TPMCOUNTS_RSEM_GENE.out.detectedgene_barplot
    ch_detectedgene_gstat_rsem_gene = CALC_TPMCOUNTS_RSEM_GENE.out.detectedgene_gstat

    //
    // MODULE: Calc detected transcripts from RSEM TPM
    //
    CALC_TPMCOUNTS_RSEM_TS (
        ch_rsem_merged_ts,
        ch_num_of_ts_rsem_header,
        ch_num_of_ts_rsem_header_gstat,
        [],
        [],
        []
    )
    ch_detectedgene_barplot_rsem_ts = CALC_TPMCOUNTS_RSEM_TS.out.detectedgene_barplot
    ch_detectedgene_gstat_rsem_ts = CALC_TPMCOUNTS_RSEM_TS.out.detectedgene_gstat

    // =================================================
    // optional: ERCC / SIRV quantification
    // =================================================

    def ch_entropy_sirv_barplot  =  Channel.empty()
    def ch_entropy_sirv_gstat  =  Channel.empty()
    if (params.spike_in_sirv) {
        //
        // MODULE: Alignment with Hisat2 [sirv]
        //
        HISAT2_SIRV (
            ch_trimmed_reads,
            ch_hisat2_sirv_idx.collect(),
            ch_tools_dir
        ).hisat2_bam_qc
        .set { ch_hisat2_bam_qc_sirv }
    
        //
        // FUNCTION: Filtering low mapped reads bams [for SIRV readcoverage]
        //
        ch_hisat2_bam_qc_sirv
        .transpose()
        .filter { name, bam, bai, flagstat -> check_mappedread_sirv(bam,flagstat,params.min_mapped_reads) }
        .map { it[0..2] }
        .set{
            ch_hisat2_bam_qc_sirv_filtered
        }
        
        //
        // MODULE: readcoverage.jl [sirv]
        //
        READCOVERAGE_SIRV (
            ch_hisat2_bam_qc_sirv_filtered
        ).readcov_sirv_results
        .set { ch_readcov_sirv_results}

        //
        // MODULE: merge sirv readcoverage outputs
        //
        MERGE_READCOVERAGE_SIRV (
            ch_readcov_sirv_results.collect()
        )

        //
        // MODULE: Quantification with RSEM [sirv]
        //
        RSEM_BOWTIE2_SIRV (
            ch_hisat2_bam_qc_sirv_filtered,
            ch_trimmed_reads,
            ch_rsem_sirv_idx.collect()
        ).rsem_isoforms_to_merge
        .set { ch_rsem_isoforms_sirv_to_merge }

        //
        // MODULE: Merge RSEM output (sirv isoforms)
        //
        MERGE_RSEM_ISOFORMS_SIRV (
            ch_rsem_isoforms_sirv_to_merge.collect()
        )
        .merged_counts
        .set{
            ch_rsem_merged_ts_sirv
        }

        //
        // MODULE: Merge RSEM output (sirv isoforms)
        //
        CALC_ENTROPY_SIRV (
            ch_rsem_merged_ts_sirv.collect(),
            ch_entropy_of_sirv_header,
            ch_entropy_of_sirv_header_gstat
        )
        ch_entropy_sirv_barplot = CALC_ENTROPY_SIRV.out.entropy_barplot
        ch_entropy_sirv_gstat = CALC_ENTROPY_SIRV.out.entropy_gstat
    }

    //=================================================
    // MODULE: MultiQC report
    //=================================================
    MULTIQC (
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_fastqc_results_raw.collect().ifEmpty([]),
        ch_fastqc_results_trimmed.collect().ifEmpty([]),
        ch_hisat2_summary.collect().ifEmpty([]),
        ch_hisat2_summary_rrna.collect().ifEmpty([]),
        ch_rseqc_results_merge.collect().ifEmpty([]),
        ch_counts_biotype.collect().ifEmpty([]),
        ch_rsem_results_stat.collect().ifEmpty([]),
        ch_sample_correlation.collect().ifEmpty([]),
        ch_ercc_correlation_barplot.collect().ifEmpty([]),
        ch_ercc_correlation_gstat.collect().ifEmpty([]),
        ch_counts_summary_all.collect().ifEmpty([]),
        ch_counts_summary_mt.collect().ifEmpty([]),
        ch_counts_summary_histone.collect().ifEmpty([]),
        ch_assignedgenome_rate_barplot.collect().ifEmpty([]),
        ch_assignedgenome_rate_gstat.collect().ifEmpty([]),
        ch_fcounts_maprate_barplot_allgene.collect().ifEmpty([]),
        ch_fcounts_maprate_gstat_allgene.collect().ifEmpty([]),
        ch_fcounts_maprate_barplot_mt.collect().ifEmpty([]),
        ch_fcounts_maprate_gstat_mt.collect().ifEmpty([]),
        ch_fcounts_maprate_barplot_histone.collect().ifEmpty([]),
        ch_fcounts_maprate_gstat_histone.collect().ifEmpty([]),
        ch_drplot.collect().ifEmpty([]),
        ch_detectedgene_barplot.collect().ifEmpty([]),
        ch_detectedgene_gstat.collect().ifEmpty([]),
        ch_detectedgene_barplot_rsem_gene.collect().ifEmpty([]),
        ch_detectedgene_gstat_rsem_gene.collect().ifEmpty([]),
        ch_detectedgene_barplot_rsem_ts.collect().ifEmpty([]),
        ch_detectedgene_gstat_rsem_ts.collect().ifEmpty([]),
        ch_entropy_sirv_barplot.collect().ifEmpty([]),
        ch_entropy_sirv_gstat.collect().ifEmpty([]),
        ch_software_versions_yaml.collect(),
        ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
    ).multiqc_report
    .set { ch_multiqc_report }
}

///////////////////////////////////////////////////////////////////////////////
/*
* workflow.onComplete
*/
///////////////////////////////////////////////////////////////////////////////

workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[ramdaq] Successful: " + workflow.runName
    if (!workflow.success) {
        subject = "[ramdaq] FAILED: " + workflow.runName
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[ramdaq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[ramdaq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[ramdaq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[ramdaq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    log.info "${c_green}----------------------------------------------------${c_reset}"
    log.info "${c_green}----------------------------------------------------${c_reset}"

    if (workflow.success) {
        log.info "${c_purple}[ramdaq]${c_green} Pipeline completed successfully!${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[ramdaq]${c_red} Pipeline completed with errors!${c_reset}"
    }
    
    // copy .nextflow.log
    today = new Date().format("yyyy-MM-dd-HH-mm-ss")
    new File("${params.outdir}/ramdaq-${today}.log") << new File('.nextflow.log').text
    
    println "${c_purple}[ramdaq]${c_green} The log file .nextflow.log was copied to ${params.outdir}/ramdaq-${today}.log${c_reset}"
}