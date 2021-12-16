// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALC_DETECTEDGENES_DR {

    label 'process_medium'
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true
    
    input:
    file(tpm_count)
    file(detectplot_header)
    file(detectplot_header_gstat)
    file(pcaplot_header)
    file(tsneplot_header)
    file(umapplot_header)
    
    output:
    path "*.{txt,pdf}", emit: detectedgene_dr_results
    path "DRplot_*.csv", optional:true, emit: drplot
    path "barplot_*.csv", emit: detectedgene_barplot
    path "gstat_*.csv", emit: detectedgene_gstat

    script:
    def prefix = "num_of_${options.suffix}"
    """
    drawplot_tpm_counts.r $tpm_count ${options.args}
    cp barplot_${prefix}.csv gstat_${prefix}.csv
    cat $detectplot_header barplot_${prefix}.csv >> tmp_file
    mv tmp_file barplot_${prefix}_mqc.csv
    cat $detectplot_header_gstat gstat_${prefix}.csv >> tmp_file
    mv tmp_file gstat_${prefix}_mqc.csv
    
    if [[ -f DRplot_pca_allsample.csv ]]; then
        cat $pcaplot_header DRplot_pca_allsample.csv >> tmp_file
        mv tmp_file DRplot_pca_allsample_mqc.csv
    fi
    
    if [[ -f DRplot_tsne_allsample.csv ]]; then
        cat $tsneplot_header DRplot_tsne_allsample.csv >> tmp_file
        mv tmp_file DRplot_tsne_allsample_mqc.csv
    fi

    if [[ -f DRplot_umap_allsample.csv ]]; then
        cat $umapplot_header DRplot_umap_allsample.csv >> tmp_file
        mv tmp_file DRplot_umap_allsample_mqc.csv
    fi
    """
}