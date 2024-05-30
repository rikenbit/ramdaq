#!/usr/bin/env nextflow
/*
========================================================================================
    ramdaq
========================================================================================
    Github : https://github.com/rikenbit/ramdaq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

/*
params.fasta         = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf           = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff           = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed      = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.bbsplit_index = WorkflowMain.getGenomeAttribute(params, 'bbsplit')
params.star_index    = WorkflowMain.getGenomeAttribute(params, 'star')
params.hisat2_index  = WorkflowMain.getGenomeAttribute(params, 'hisat2')
params.rsem_index    = WorkflowMain.getGenomeAttribute(params, 'rsem')
params.salmon_index  = WorkflowMain.getGenomeAttribute(params, 'salmon')
*/
      /*
      adapter = "${params.local_annot_dir}/all_sequencing_WTA_adopters.fa"
      hisat2_idx = "${params.local_annot_dir}/hisat2_v220_index_GRCm39_primary_ERCC/*.ht2"
      hisat2_rrna_idx = "${params.local_annot_dir}/hisat2_index_Mouse_rRNA_12_mito_5_45.tar.gz"
      chrsize = "${params.local_annot_dir}/GRCm39.primary_assembly.genome.fa.fai"
      bed = "${params.local_annot_dir}/gencode.vM26.primary_assembly.annotation.bed"
      gtf = "${params.local_annot_dir}/gencode.vM26.primary_assembly.annotation.ERCC.gtf"
      mt_gtf = "${params.local_annot_dir}/gencode.vM26.primary_assembly.annotation.mt.gtf"
      rrna_gtf = "${params.local_annot_dir}/Mouse_GRCm39_rmsk_rRNA_Refseq_45S_merge.gtf"
      histone_gtf = "${params.local_annot_dir}/gencode.vM26.primary_assembly.annotation.histone.gtf"
      hisat2_sirv_idx = "${params.local_annot_dir}/hisat2_v220_index_SIRVome/*.ht2"
      rsem_sirv_idx = "${params.local_annot_dir}/RSEM_bowtie2_index_SIRVome/*"
      rsem_allgene_idx = "${params.local_annot_dir}/RSEM_bowtie2_index_gencode_vM26/*"
      */

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

//WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { DL_REFERENCES } from './workflows/download_references'
include { RAMDAQ }        from './workflows/ramdaq'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    
    if (params.dl_references) {
        DL_REFERENCES()
    } else {
        RAMDAQ()
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/