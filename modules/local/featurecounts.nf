// Import generic module functions
include { initOptions } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FEATURECOUNTS {
    tag "$name"
    
    publishDir "${params.outdir}/${options.publish_dir}", mode: 'copy', overwrite: true,
        saveAs: {filename ->
            if (filename.indexOf("biotype_counts") > 0) "biotype_counts/$filename"
            else if (filename.indexOf(".featureCounts.txt.summary") > 0) "count_summaries/$filename"
            else if (filename.indexOf(".featureCounts.txt") > 0) "counts/$filename"
            else "$filename"
        }
    
    input:
    tuple val(name), file(bam), file(bai)
    file gtf
    file biotypes_header

    output:
    path "*.featureCounts.txt", emit: counts_to_merge
    path "*.featureCounts.txt.summary", emit: counts_summary
    path "${name}.allgene.featureCounts.txt", optional:true, emit: counts_to_plot_corr
    path "${name}_biotype_counts*mqc.{txt,tsv}", optional:true, emit: counts_biotype

    script:
    def prefix = options.suffix ? "${name}${options.suffix}" : "${name}"

    def is_pairedend = params.single_end ? '' : "-p"
    def strandspecific = ''
    if (params.stranded && params.stranded == 'fr-firststrand') {
        strandspecific = "-s 2"
    } else if (params.stranded && params.stranded == 'fr-secondstrand'){
        strandspecific = "-s 1"
    }
    def extra_attributes = params.extra_attributes ? "--extraAttributes ${params.extra_attributes}" : ''
    def allow_multimap = params.allow_multimap ? "-M" : ''
    def allow_overlap = params.allow_overlap ? "-O" : ''
    def count_fractionally = params.count_fractionally ? "--fraction" : ''
    def threads_num = params.fc_threads_num > 0 ? "-T ${params.fc_threads_num}" : ''
    def biotype = params.group_features_type

    if (options.suffix != '.allgene') {
        """
        featureCounts -a $gtf -g ${params.group_features} -t ${params.count_type} -o ${prefix}.featureCounts.txt  \\
        $is_pairedend $strandspecific $extra_attributes $allow_multimap $allow_overlap $count_fractionally $threads_num ${bam}
        """
    } else {
        """
        featureCounts -a $gtf -g ${params.group_features} -t ${params.count_type} -o ${prefix}.featureCounts.txt  \\
        $is_pairedend $strandspecific $extra_attributes $allow_multimap $allow_overlap $count_fractionally $threads_num ${bam}

        featureCounts -a $gtf -g $biotype -o ${name}_biotype_featureCounts.txt $is_pairedend $strandspecific ${bam}
        cut -f 1,7 ${name}_biotype_featureCounts.txt | tail -n +3 | cat $biotypes_header - >> ${name}_biotype_counts_mqc.txt
        """

    }




}