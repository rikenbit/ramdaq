#!/usr/bin/env Rscript

# params : NuclearRNA genes
nuclearRNA_human = c("MALAT1", "NEAT1", "RN7SK", "RN7SL1")
nuclearRNA_mouse = c("Malat1", "Neat1", "Rn7sk", "Rn7s1")

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

if (!require("stringr")){
    install.packages("stringr", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("stringr")
}
if (!require("dplyr")){
    install.packages("dplyr", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("dplyr")
}
if (!require("ggplot2")){
    install.packages("ggplot2", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("ggplot2")
}

inputfile01 <- args[1]
raw_countdata = read.table(inputfile01, sep="\t", comment.char = "", header=T, check.names=FALSE, stringsAsFactors=F)

inputfile02 <- args[2]
tpm_countdata = read.table(inputfile02, sep="\t", comment.char = "", header=T, check.names=FALSE, stringsAsFactors=F)

nuclearRNA_qc <- args[3]

nuclearRNA_genes = c()
if (nuclearRNA_qc == "mouse") {
  nuclearRNA_genes = nuclearRNA_mouse
} else if (nuclearRNA_qc == "human"){
  nuclearRNA_genes = nuclearRNA_human
}

# create gene data
gene_data = data.frame(Geneid=raw_countdata$Geneid, gene_name=raw_countdata$gene_name, stringsAsFactors=F)

# create tpm data
rownames(tpm_countdata) = tpm_countdata[,1]
if (ncol(tpm_countdata) > 3){
  tpm_countdata = tpm_countdata[,-c(1,2)]
} else {
  tmp_rowname = rownames(tpm_countdata)
  tmp_colname = colnames(tpm_countdata)
  tpm_countdata = data.frame(tmp = tpm_countdata[,-c(1,2)])
  colnames(tpm_countdata) = tmp_colname[3]
  rownames(tpm_countdata) = tmp_rowname
}

# create nuclear-RNA counts
if (length(nuclearRNA_genes) > 0) {

    nuclearRNA_geneid = subset(raw_countdata, gene_name %in% nuclearRNA_genes)$Geneid

    if (ncol(tpm_countdata) > 1){
      count_tpm_nuclearRNA = tpm_countdata[rownames(tpm_countdata) %in% nuclearRNA_geneid,]
    } else {
      count_tpm_nuclearRNA = data.frame(tmp = tpm_countdata[rownames(tpm_countdata) %in% nuclearRNA_geneid,])
      colnames(count_tpm_nuclearRNA) = colnames(tpm_countdata)
      rownames(count_tpm_nuclearRNA) = nuclearRNA_geneid
    }
    
    count_tpm_nuclearRNA_log = log10(count_tpm_nuclearRNA+1)
    count_tpm_nuclearRNA_log$Geneid = rownames(count_tpm_nuclearRNA_log)
    count_tpm_nuclearRNA_log = dplyr::left_join(count_tpm_nuclearRNA_log, gene_data, by=c("Geneid"))
    write.table(count_tpm_nuclearRNA_log,"merged_featureCounts_nuclearRNA_log.txt",sep="\t", append=F, quote=F, col.names=T)

    outputdata = count_tpm_nuclearRNA_log[,!colnames(count_tpm_nuclearRNA_log) %in% c("Geneid", "gene_name")]
    outputdata = as.data.frame(t(outputdata))
    colnames(outputdata) = count_tpm_nuclearRNA_log$gene_name
    outputdata$name = colnames(count_tpm_nuclearRNA_log)[!colnames(count_tpm_nuclearRNA_log) %in% c("Geneid", "gene_name")]
    rownames(outputdata) = colnames(count_tpm_nuclearRNA_log)[!colnames(count_tpm_nuclearRNA_log) %in% c("Geneid", "gene_name")]
    outputdata = outputdata[,c("name", count_tpm_nuclearRNA_log$gene_name)]
    write.csv(outputdata, 'barplot_nuclear_rna_exp.csv', quote=FALSE, append=TRUE)
}

# Printing sessioninfo to standard out
print("Draw plot info:")
sessionInfo()







