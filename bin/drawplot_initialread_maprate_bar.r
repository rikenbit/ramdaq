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

input_seqcounts <- args[1]
input_totalread <- args[2]
isPairedEnd <- args[3]

### function ###
seqcountsdata = read.table(input_seqcounts, sep="", comment.char = "", header=F, check.names=FALSE, stringsAsFactors=F)
totalreaddata = read.table(input_totalread, sep="", comment.char = "", header=F, check.names=FALSE, stringsAsFactors=F)

plotdata = dplyr::left_join(seqcountsdata, totalreaddata, by=c("V1"))
colnames(plotdata) = c("samplename", "initialreads", "totalread")
plotdata[is.na(plotdata)] <- 0

if (isPairedEnd=="True"){
    plotdata$initialReadsMapRate = (plotdata$totalread / (plotdata$initialreads*2)) *100
} else {
    plotdata$initialReadsMapRate = (plotdata$totalread / plotdata$initialreads) *100
}

plotdata = plotdata[order(plotdata$samplename),]

### draw barplot
g = ggplot(plotdata, aes(x=samplename,y=initialReadsMapRate)) +
    geom_bar(alpha=0.7, stat="identity") +
    ylim(0, 100) +
    xlab("Sample")+ylab("Initial reads mapped rate") +
    #scale_x_discrete(labels=plotdata$x) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle("Rate of initial reads mapped to genome")

ggsave(file = "barplot_initialread_maprate.pdf", plot=g, dpi=100, width=12, height=5)

# Write rate values to file
outputdata = data.frame(name = plotdata$samplename, initialReadsMapRate=plotdata$initialReadsMapRate)
rownames(outputdata) = outputdata$name
write.csv(outputdata, 'barplot_initialread_maprate.csv', quote=FALSE, append=TRUE)

# Printing sessioninfo to standard out
print("Draw plot info:")
sessionInfo()

