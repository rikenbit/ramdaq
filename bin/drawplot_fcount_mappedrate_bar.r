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
if (!require("ggplot2")){
    install.packages("ggplot2", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("ggplot2")
}

input_totalseq <- args[1]
input_fcount_merged <- args[2]
isPairedEnd <- args[3]
annotation_name <- args[4]

### function ###
totalseq_data = read.table(input_totalseq, sep="", comment.char = "", header=F, check.names=FALSE, stringsAsFactors=F)
fcount_merged_data = read.table(input_fcount_merged, sep="", comment.char = "", header=T, check.names=FALSE, stringsAsFactors=F)

if (ncol(fcount_merged_data) > 4){
  fcount_merged_data = fcount_merged_data[,!colnames(fcount_merged_data) %in% c("Geneid", "Length", "gene_name")]
} else {
  tmp_rowname = fcount_merged_data$Geneid
  tmp_colname = colnames(fcount_merged_data)
  fcount_merged_data = data.frame(tmp = fcount_merged_data[,-c(1,2,3)])
  colnames(fcount_merged_data) = tmp_colname[4]
  rownames(fcount_merged_data) = tmp_rowname
}

total_mapped_reads = data.frame(V1 = colnames(fcount_merged_data), mappedreads = colSums(fcount_merged_data), stringsAsFactors=F)

plotdata = dplyr::left_join(totalseq_data, total_mapped_reads, by=c("V1"))
colnames(plotdata) = c("samplename", "totalseq", "mappedreads")
plotdata[is.na(plotdata)] <- 0

if (isPairedEnd=="True"){
    plotdata$assignedRate = (plotdata$mappedreads / (plotdata$totalseq*2)) *100
} else {
    plotdata$assignedRate = (plotdata$mappedreads / plotdata$totalseq) *100
}

plotdata = plotdata[order(plotdata$samplename),]

### draw barplot
g = ggplot(plotdata, aes(x=samplename,y=assignedRate)) +
    geom_bar(alpha=0.7, stat="identity") +
    ylim(0, 100) +
    xlab("sample")+ylab("assigned rate") +
    #scale_x_discrete(labels=plotdata$x) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle(paste0("Rate of reads assigned to ", annotation_name, " annotation"))

ggsave(file = paste0("barplot_assigned_",annotation_name,"_rate.pdf"), plot=g, dpi=100, width=12, height=5)

# Write rate values to file
outputdata = data.frame(name = plotdata$samplename, assignedRate=plotdata$assignedRate)
rownames(outputdata) = outputdata$name
write.csv(outputdata, paste0("barplot_assignedrate_",annotation_name,".csv"), quote=FALSE, append=TRUE)

# Printing sessioninfo to standard out
print("Draw plot info:")
sessionInfo()



