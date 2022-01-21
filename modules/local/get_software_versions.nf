// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_SOFTWARE_VERSIONS {

    label 'process_low'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true,
        saveAs: { filename ->
                  if (filename.indexOf(".csv") > 0) filename
                  else null
            }
    
    output:
    path "software_versions_mqc.yaml", emit: software_versions_yaml
    path "software_versions.csv"
    
    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    fastq-mcf -V > v_fastqmcf.txt
    hisat2 --version > v_hisat2.txt
    samtools --version > v_samtools.txt
    bam2wig.py --version > v_bam2wig.txt
    bamtools --version > v_bamtools.txt
    read_distribution.py --version > v_read_distribution.txt
    infer_experiment.py --version > v_infer_experiment.txt
    inner_distance.py --version > v_inner_distance.txt
    junction_annotation.py --version > v_junction_annotation.txt
    featureCounts -v > v_featurecounts.txt 2>&1
    rsem-calculate-expression --version > v_rsem.txt
    Rscript -e "write(x=as.character(R.version.string), file='v_R.txt')"
    Rscript -e "library(edgeR); write(x=as.character(packageVersion('edgeR')), file='v_edgeR.txt')"
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}