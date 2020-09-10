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

ercc_countdata = read.table(erccfile, sep="\t", comment.char = "", header=T, stringsAsFactors=F)
ercc_ref = read.table(erccref, sep="\t", comment.char = "", header=T, stringsAsFactors=F)

plotdata = cor(ercc_countdata, ercc_ref$mollog)
plotdata = data.frame(corr=plotdata[,1],sample_name=as.character(rownames(plotdata)))
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





