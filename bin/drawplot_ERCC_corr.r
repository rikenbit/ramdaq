#!/usr/bin/env Rscript

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

if (!require("ggplot2")){
    install.packages("ggplot2", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("ggplot2")
}

erccfile <- args[1]
erccref <- args[2]

### function ###

SDset <- function(x, y){
  data.frame(x=x, y=y, upper = mean(y)+3*sd(y),lower = mean(y)-3*sd(y))
}

calc_tpm <- function(counts,len) {
  x <- counts/len
    return(t(t(x)*1e6/colSums(x)))
}

### load Ref
ercc_ref = read.table(erccref, sep="\t", comment.char = "", header=T, stringsAsFactors=F)
ercc_ref_trim = ercc_ref[,c("ERCC.ID","Length")]

ercc_countdata = read.table(erccfile, sep="\t", comment.char = "", header=T, stringsAsFactors=F)
ercc_countdata$ERCC.ID = rownames(ercc_countdata)
ercc_countdata_tpm = dplyr::left_join(ercc_countdata, ercc_ref_trim, by=c("ERCC.ID"))
rownames(ercc_countdata_tpm) = ercc_countdata_tpm$ERCC.ID

# calc TPM
ercc_countdata_tpm = calc_tpm(ercc_countdata_tpm[,!colnames(ercc_countdata_tpm) %in% c("ERCC.ID","Length"), drop=FALSE], as.numeric(ercc_countdata_tpm$Length))
ercc_countdata_tpm = as.data.frame(ercc_countdata_tpm)
ercc_countdata_tpm_log = log10(ercc_countdata_tpm+1)

if (nrow(ercc_countdata_tpm_log) > 1){
  write.table(ercc_countdata_tpm_log,"merged_featureCounts_gene_TPM_ERCC_log.txt",sep="\t", append=F, quote=F, row.names=T, col.names=T)
}

plotdata = cor(ercc_countdata_tpm_log, ercc_ref$mollog)
plotdata = data.frame(corr=plotdata[,1],sample_name=as.character(rownames(plotdata)))
plotdata = plotdata[order(plotdata$sample_name),]
plotdata = SDset(plotdata$sample_name, plotdata$corr)

### draw plot
g = ggplot(plotdata, aes(x=x,y=y)) +
    geom_bar(alpha=0.7, stat="identity") +
    xlab("Sample") + ylab("Corr") + 
    ylim(0,1) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle("ERCC counts and mol-log correlation")

ggsave(file = "ERCC_countsmol_correlation.pdf", plot=g, dpi=100, width=12, height=5)

# Write correlation values to file
output_data = data.frame(name = plotdata$x, corr=plotdata$y)
rownames(output_data) = output_data$name
write.csv(output_data, 'ercc_countsmol_correlation.csv', quote=FALSE, append=TRUE)

# Printing sessioninfo to standard out
print("Draw ERCC correlation plot info:")
sessionInfo()





