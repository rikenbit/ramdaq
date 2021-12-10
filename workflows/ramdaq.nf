
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


/*
========================================================================================
    SET UP REPORT HEADERS
========================================================================================
*/

ch_biotypes_header = file("$projectDir/assets/biotypes_header.txt", checkIfExists: true)

/*
========================================================================================
    SET UP DIR/FILE PATH VARIABLES
========================================================================================
*/

// Tools dir
ch_tools_dir = workflow.scriptFile.parent + "/tools"

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

// Function that checks the number of mapped reads from flagstat output
// and returns true if > params.min_mapped_reads and otherwise false
def check_mappedread(name,flagstat_path,min_mapped_reads=10) {
    def flagstat_file = new File(flagstat_path.toString())
    def mapped = 0
    flagstat_file.eachLine { line ->
        if (line.contains(' mapped (')) {
            mapped = line.tokenize().first().toInteger()
        }
    }
    if (mapped < min_mapped_reads.toInteger()) {
        log.info ">>>> $name FAILED MAPPED READ THRESHOLD: ${mapped} < ${params.min_mapped_reads}. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! <<<<"
        return false
    } else {
        return true
    }
}

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { UNTAR_INDEX as UNTAR_HISAT2_IDX } from '../modules/local/untar_index' addParams( options: modules['untar_index_hisat2'] )
include { UNTAR_INDEX as UNTAR_HISAT2_RRNA_IDX } from '../modules/local/untar_index' addParams( options: modules['untar_index_hisat2'] )
include { UNTAR_INDEX as UNTAR_RSEM_ALLGENE_IDX } from '../modules/local/untar_index' addParams( options: modules['untar_index_rsem'] )

include { FASTQC as FASTQC_RAW } from '../modules/local/fastqc' addParams( options: modules['fastqc'] )
include { FASTQC as FASTQC_TRIM } from '../modules/local/fastqc' addParams( options: modules['fastqc_trim'] )
include { FASTQMCF } from '../modules/local/fastqmcf' addParams( options: modules['fastqmcf'] )

include { HISAT2 as HISAT2_GENOME } from '../modules/local/hisat2' addParams( options: modules['hisat2_genome'] )
include { HISAT2 as HISAT2_RRNA } from '../modules/local/hisat2' addParams( options: modules['hisat2_rrna'] )
include { MERGE_SUMMARYFILE as MERGE_SUMMARYFILE_HISAT2 } from '../modules/local/merge_summaryfile' addParams( options: modules['merge_summaryfile_hisat2'] )

include { BAM2WIG as BAM2WIG_ALLGENES } from '../modules/local/bam2wig' addParams( options: modules['bam2wig'] )
include { ADJUST_BED_NONCODING } from '../modules/local/adjust_bed_noncoding' addParams( options: modules['adjust_bed_noncoding'] )
include { RSEQC } from '../modules/local/rseqc' addParams( options: modules['rseqc'] )
include { READCOVERAGE } from '../modules/local/readcoverage' addParams( options: modules['readcoverage'] )

include { FEATURECOUNTS as FEATURECOUNTS_ALL_GTF} from '../modules/local/featurecounts' addParams( options: modules['featurecounts_all_gtf'] )
include { FEATURECOUNTS as FEATURECOUNTS_MT_GTF} from '../modules/local/featurecounts' addParams( options: modules['featurecounts_mt_gtf'] )
include { FEATURECOUNTS as FEATURECOUNTS_HISTONE_GTF} from '../modules/local/featurecounts' addParams( options: modules['featurecounts_histone_gtf'] )
include { MERGE_FEATURECOUNTS as MERGE_FEATURECOUNTS_ALLGENE} from '../modules/local/merge_featurecounts' addParams( options: modules['merge_featurecounts_allgene'] )
include { MERGE_FEATURECOUNTS as MERGE_FEATURECOUNTS_MT} from '../modules/local/merge_featurecounts' addParams( options: modules['merge_featurecounts_mt'] )
include { MERGE_FEATURECOUNTS as MERGE_FEATURECOUNTS_HISTONE} from '../modules/local/merge_featurecounts' addParams( options: modules['merge_featurecounts_histone'] )

include { RSEM_BOWTIE2 as RSEM_BOWTIE2_ALLGENES } from '../modules/local/rsem_bowtie2' addParams( options: modules['rsem_bowtie2_allgenes'] )
include { MERGE_RSEM as MERGE_RSEM_GENES} from '../modules/local/merge_rsem' addParams( options: modules['merge_rsem_genes'] )
include { MERGE_RSEM as MERGE_RSEM_ISOFORMS} from '../modules/local/merge_rsem' addParams( options: modules['merge_rsem_isoforms'] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow RAMDAQ {
    //ch_versions = Channel.empty()

    if (!ch_hisat2_idx){
        //
        // MODULE: untar index.tar.gz [hisat2 all genes]
        //
        UNTAR_HISAT2_IDX (
            params.hisat2_idx
        )
        .index_files
        .set { ch_hisat2_idx }
    }
    if (!ch_hisat2_rrna_idx){
        //
        // MODULE: untar index.tar.gz [hisat2 rrna]
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
    // FUNCTION: Filtering low mapped reads bams [for bamQC]
    //
    ch_hisat2_bam_qc
    .transpose()
    .filter { name, bam, bai, flagstat -> check_mappedread(name,flagstat,params.min_mapped_reads) }
    .map { it[0..2] }
    .set{
        ch_hisat2_bam_qc_filtered
    }

    //
    // FUNCTION: Filtering low mapped reads bams [for featureCounts]
    //
    ch_hisat2_bam_count
    .filter { name, bam, bai, flagstat -> check_mappedread(name,flagstat,params.min_mapped_reads) }
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
    .filter { name, bam, bai, flagstat -> check_mappedread(bam.baseName,flagstat,params.min_mapped_reads) }
    .map { it[1] }
    .set { ch_bam_samplenum }

    //
    // MODULE: Merge summaryfiles [Hisat2 totalseq]
    //
    MERGE_SUMMARYFILE_HISAT2 (
        ch_hisat2_summary.collect()
    )

    //
    // MODULE: Alignment with Hisat2 [rrna]
    //
    HISAT2_RRNA (
        ch_trimmed_reads,
        ch_hisat2_rrna_idx.collect(),
        ch_tools_dir
    )

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
    .readdist_totalread
    .set { ch_readdist_totalread }

    //
    // MODULE: readcoverage.jl
    //
    READCOVERAGE (
        ch_hisat2_bam_qc_filtered,
        ch_bed
    )

    //
    // MODULE: featureCounts (All-genes GTF)
    //
    FEATURECOUNTS_ALL_GTF (
        ch_hisat2_bam_featurecount,
        ch_gtf,
        ch_biotypes_header
    )
    ch_counts_to_merge_all = FEATURECOUNTS_ALL_GTF.out.counts_to_merge
    ch_counts_summary_all = FEATURECOUNTS_ALL_GTF.out.counts_summary
    ch_counts_to_plot_corr = FEATURECOUNTS_ALL_GTF.out.counts_to_plot_corr

    //
    // MODULE: Merge featureCounts output (All-genes)
    //
    MERGE_FEATURECOUNTS_ALLGENE (
        ch_counts_to_merge_all.collect()
    )

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

    //
    // MODULE: Quantification with RSEM [all genes]
    //
    RSEM_BOWTIE2_ALLGENES (
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

    //
    // MODULE: Merge RSEM output (isoforms)
    //
    MERGE_RSEM_ISOFORMS (
        ch_rsem_isoforms_to_merge.collect()
    )
}

