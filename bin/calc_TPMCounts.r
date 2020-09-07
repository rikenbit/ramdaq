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
raw_countdata = read.table(inputfile, sep="\t", comment.char = "", header=T, stringsAsFactors=F)

# create gene length data
gene_length = data.frame(Geneid=raw_countdata$Geneid, Length=raw_countdata$Length, stringsAsFactors=F)

### function ###

trim_countfc_sample <- function(data, ERCC=F){

  if (ERCC) {
    data = subset(data, grepl("ERCC", Geneid))
  } else {
    data = subset(data, !grepl("ERCC", Geneid))
  }

  rownames(data) = data$Geneid
  data = data[,-c(1,2,3), drop=FALSE]
  samplename = str_replace(colnames(data), ".sort", "")
  colnames(data) = c(samplename)
  return(data)
}

calc_tpm <- function(counts,len) {
  x <- counts/len
    return(t(t(x)*1e6/colSums(x)))
}

count_trim = trim_countfc_sample(raw_countdata)
count_trim$Geneid = rownames(count_trim)
count_tpm = dplyr::left_join(count_trim, gene_length, by=c("Geneid"))
rownames(count_tpm) = count_tpm$Geneid

# calc TPM
count_tpm = calc_tpm(count_tpm[,!colnames(count_tpm) %in% c("Geneid","Length"), drop=FALSE], as.numeric(count_tpm$Length))
count_tpm = as.data.frame(count_tpm)

# Print distance matrix to file
write.table(count_tpm,"merged_featureCounts_gene_TPM.txt",sep="\t", append=F, quote=F, row.names=T, col.names=T)

# create ERCC counts
count_ercc = trim_countfc_sample(raw_countdata, ERCC=T)

if (nrow(count_ercc) > 1){
  write.table(count_ercc,"merged_featureCounts_gene_ERCC.txt",sep="\t", append=F, quote=F, row.names=T, col.names=T)
}

# Printing sessioninfo to standard out
print("Calc TPM counts info:")
sessionInfo()





