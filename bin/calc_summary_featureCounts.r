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
merged_data = read.table(inputfile, sep="\t", comment.char = "", header=T, check.names=FALSE, stringsAsFactors=F)
merged_data_t = t(merged_data[,-1])
colnames(merged_data_t) = merged_data[,1]
merged_data_t = as.data.frame(merged_data_t)
merged_data_t$Total = rowSums(merged_data_t)
merged_data_t$percent_assigned = merged_data_t$Assigned / merged_data_t$Total

write.table(merged_data_t,"merged_featureCounts_summary_rrna.txt",sep="\t", append=F, quote=F, row.names=T, col.names=T)

# Printing sessioninfo to standard out
print("Calc featureCounts summary info:")
sessionInfo()





