#!/usr/bin/env Rscript

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

inputfile <- args[1]
raw_countdata = read.table(inputfile, sep="\t", comment.char = "", header=T, check.names=FALSE, stringsAsFactors=F)

### function ###

trim_tpmmatrix <- function(data, ERCC=F){

  if (ncol(data) > 1){
    if (ERCC) {
      data = data[grepl("ERCC", rownames(data)),]
    } else {
      data = data[!grepl("ERCC", rownames(data)),]
    }
  } else {
    tmp_rowname = rownames(data)
    tmp_colname = colnames(data)
    if (ERCC) {
      data = data.frame( tmp = data[grepl("ERCC", rownames(data)),])
      tmp_rowname = tmp_rowname[grepl("ERCC", tmp_rowname)]
    } else {
      data = data.frame( tmp = data[!grepl("ERCC", rownames(data)),])
      tmp_rowname = tmp_rowname[!grepl("ERCC", tmp_rowname)]
    }
    colnames(data) = tmp_colname
    rownames(data) = tmp_rowname
  }

  return(data)
}

calc_tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

count_tpm = raw_countdata

# calc TPM
count_tpm = calc_tpm(count_tpm[,!colnames(count_tpm) %in% c("Geneid","Length","gene_name"), drop=FALSE], as.numeric(count_tpm$Length))
count_tpm = as.data.frame(count_tpm)
rownames(count_tpm) = raw_countdata$Geneid

# replace NaN value with zero
count_tpm = as.matrix(count_tpm)
count_tpm[is.nan(count_tpm)] <- 0
count_tpm = as.data.frame(count_tpm)

# create all-genes counts
count_tpm_gene = trim_tpmmatrix(count_tpm)
write.table(count_tpm_gene,"merged_featureCounts_allgene_TPM.txt",sep="\t", append=F, quote=F, row.names=T, col.names=T)

# create ERCC counts
count_tpm_ercc = trim_tpmmatrix(count_tpm, ERCC=T)

if (nrow(count_tpm_ercc) > 1 && sum(colSums(count_tpm_ercc)) > 0){
  count_tpm_ercc_log = log10(count_tpm_ercc+1)
  write.table(count_tpm_ercc_log,"merged_featureCounts_ERCC_TPM_log.txt",sep="\t", append=F, quote=F, row.names=T, col.names=T)

  raw_countdata_ercc = subset(raw_countdata, grepl("ERCC", Geneid))
  rownames(raw_countdata_ercc) = raw_countdata_ercc$Geneid
  raw_countdata_ercc = raw_countdata_ercc[,!colnames(raw_countdata_ercc) %in% c("Geneid","Length","gene_name"), drop=FALSE]
  write.table(raw_countdata_ercc,"merged_featureCounts_ERCC.txt",sep="\t", append=F, quote=F, row.names=T, col.names=T)
}

# Printing sessioninfo to standard out
print("Calc TPM counts info:")
sessionInfo()





