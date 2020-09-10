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
input_totalread <- args[2]
isPairedEnd <- args[3]

### function ###
totalseqdata = read.table(input_totalseq, sep="", comment.char = "", header=F, stringsAsFactors=F)
totalreaddata = read.table(input_totalread, sep="", comment.char = "", header=F, stringsAsFactors=F)

plotdata = dplyr::left_join(totalseqdata, totalreaddata, by=c("V1"))
colnames(plotdata) = c("samplename", "totalseq", "totalread")

if (isPairedEnd=="True"){
    plotdata$assignedGeneRate = (plotdata$totalread / (plotdata$totalseq*2)) *100
} else {
    plotdata$assignedGeneRate = (plotdata$totalread / plotdata$totalseq) *100
}

### draw barplot
g = ggplot(plotdata, aes(x=samplename,y=assignedGeneRate)) +
    geom_bar(alpha=0.7, stat="identity") +
    ylim(0, 100) +
    xlab("Sample")+ylab("Assigned gene rate") +
    #scale_x_discrete(labels=plotdata$x) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle("Rate of assigned to genes")

ggsave(file = "barplot_assignedgene_rate.pdf", plot=g, dpi=100, width=12, height=5)

# Write rate values to file
outputdata = data.frame(name = plotdata$samplename, assignedGeneRate=plotdata$assignedGeneRate)
rownames(outputdata) = outputdata$name
write.csv(outputdata, 'barplot_assignedgene_rate.csv', quote=FALSE, append=TRUE)

# Printing sessioninfo to standard out
print("Draw plot info:")
sessionInfo()







