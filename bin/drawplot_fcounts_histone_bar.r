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

input_histone <- args[1]

### function ###
histone_summary= read.table(input_histone, sep="", comment.char = "", header=T, stringsAsFactors=F)
rownames(histone_summary) = histone_summary[,1]
histone_summary = as.data.frame(t(histone_summary[,-1]))
histone_summary$Total_Counts = rowSums(histone_summary)
histone_summary$Assigned_Rate =  (histone_summary$Assigned/histone_summary$Total_Counts)*100
histone_summary$Sample_Name = rownames(histone_summary)

plotdata = histone_summary[,c("Sample_Name","Assigned_Rate")]
colnames(plotdata) = c("samplename", "assignedHistoneRate")
plotdata = plotdata[order(plotdata$samplename),]

### draw barplot
g = ggplot(plotdata, aes(x=samplename,y=assignedHistoneRate)) +
    geom_bar(alpha=0.7, stat="identity") +
    #ylim(0, 100) +
    xlab("Sample")+ylab("Assigned Histone genes rate") +
    #scale_x_discrete(labels=plotdata$x) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle("Rate of assigned to Histone genes")

ggsave(file = "barplot_assigned_histone_rate.pdf", plot=g, dpi=100, width=12, height=5)

# Write rate values to file
outputdata = data.frame(name = plotdata$samplename, assignedHistoneRate=plotdata$assignedHistoneRate)
rownames(outputdata) = outputdata$name
write.csv(outputdata, 'barplot_fcounts_histone_rate.csv', quote=FALSE, append=TRUE)

# Printing sessioninfo to standard out
print("Draw plot info:")
sessionInfo()







