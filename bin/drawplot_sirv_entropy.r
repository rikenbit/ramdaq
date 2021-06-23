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

### function ###
shannon.entropy <- function(p)
{
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

input_tpm <- args[1]

### create sirv counts
sirv_tpm = read.table(input_tpm, sep="", comment.char = "", header=T, check.names=FALSE, stringsAsFactors=F)

colnames(sirv_tpm) = str_split_fixed(colnames(sirv_tpm), ".sirv", 2)[,1]
sirv_tpm = as.data.frame(sirv_tpm)
sirv_tpm = subset(sirv_tpm, !grepl("ERCC", gene_id))

rownames(sirv_tpm) = sirv_tpm$transcript_id
sirv_tpm = sirv_tpm[,-c(1,2)]

### calc entropy
samplenum = ncol(sirv_tpm)
sirvnum = nrow(sirv_tpm)
sirv_entropies = c()
for(j in 1:samplenum){
  p <- sirv_tpm[1:sirvnum,j]/sum(sirv_tpm[1:sirvnum,j])
  sirv_entropies[j] <- ifelse(!any(is.na(p)), shannon.entropy(p), -1)
}

plotdata = data.frame(samplename = colnames(sirv_tpm), entropy = sirv_entropies)

### draw barplot

ymax = shannon.entropy(rep(1/sirvnum, sirvnum))

g = ggplot(plotdata, aes(x=samplename,y=entropy)) +
    geom_bar(alpha=0.7, stat="identity") +
    ylim(0, ymax) +
    xlab("Sample")+ylab("Entropy of sirv") +
    #scale_x_discrete(labels=plotdata$x) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle("Entropy of sirv isoforms")

ggsave(file = "barplot_entropy_of_sirv.pdf", plot=g, dpi=100, width=12, height=5)

# Write rate values to file
outputdata = data.frame(name = plotdata$samplename, entropy=plotdata$entropy)
rownames(outputdata) = outputdata$name
write.csv(outputdata, 'barplot_entropy_of_sirv.csv', quote=FALSE, append=TRUE)

# Printing sessioninfo to standard out
print("Draw plot info:")
sessionInfo()


