#!/usr/bin/env Rscript

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

if (!require("ggplot2")){
    install.packages("ggplot2", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("ggplot2")
}

erccfile <- args[1]
erccref <- args[2]
erccamount <- args[3]

### function ###

SDset <- function(x, y){
  data.frame(x=x, y=y, upper = mean(y)+3*sd(y),lower = mean(y)-3*sd(y))
}

calc_tpm <- function(counts,len) {
  x <- counts/len
    return(t(t(x)*1e6/colSums(x)))
}

### load Ref
ercc_ref = read.table(erccref, sep="\t", comment.char = "", header=T, check.names=FALSE, stringsAsFactors=F)

### calc copy number
ercc_input_amount = eval(parse(text = as.character(erccamount)))
ercc_ref$copy.number = ercc_ref$concentration.in.Mix.1..attomoles.ul. * 10^(-18) * 6.02214086  * 10^23 * ercc_input_amount
ercc_ref$copy.number.log = log10(ercc_ref$copy.number+1)

write.table(ercc_ref,"ercc_dataset_user.txt",sep="\t", append=F, quote=F, row.names=T, col.names=T)

### load countdata
ercc_countdata_tpm_log = read.table(erccfile, sep="\t", comment.char = "", header=T, check.names=FALSE, stringsAsFactors=F)
ercc_ref_sort = ercc_countdata_tpm_log
ercc_ref_sort$ERCC.ID = rownames(ercc_ref_sort)
ercc_ref_sort = dplyr::left_join(ercc_ref_sort, ercc_ref, by=c("ERCC.ID"))

# calc corr
plotdata = cor(ercc_countdata_tpm_log, ercc_ref_sort$copy.number.log)
plotdata = data.frame(corr=plotdata[,1],sample_name=as.character(rownames(plotdata)))
plotdata = plotdata[order(plotdata$sample_name),]
plotdata = SDset(plotdata$sample_name, plotdata$corr)

### draw plot
g = ggplot(plotdata, aes(x=x,y=y)) +
    geom_bar(alpha=0.7, stat="identity") +
    xlab("Sample") + ylab("Corr") + 
    ylim(0,1) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle("ERCC counts and copy number correlation")

ggsave(file = "barplot_ercc_counts_copynum_correlation.pdf", plot=g, dpi=100, width=12, height=5)

# Write correlation values to file
output_data = data.frame(name = plotdata$x, corr=plotdata$y)
rownames(output_data) = output_data$name
write.csv(output_data, 'barplot_ercc_counts_copynum_correlation.csv', quote=FALSE, append=TRUE)

# Printing sessioninfo to standard out
print("Draw ERCC correlation plot info:")
sessionInfo()





