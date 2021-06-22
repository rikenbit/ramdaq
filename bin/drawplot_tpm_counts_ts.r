#!/usr/bin/env Rscript

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

if (!require("stringr")){
    install.packages("stringr", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("stringr")
}

if (!require("ggplot2")){
    install.packages("ggplot2", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("ggplot2")
}

inputfile <- args[1]

### function ###

SDset <- function(x, y){
  data.frame(x=x, y=y, upper = mean(y)+3*sd(y),lower = mean(y)-3*sd(y))
}

ggplot_2D <- function(dset, x, y, title, celltype, outname, label=F) {
  .e = environment()
  p=ggplot(dset, aes_string(x=x, y=y), environment=.e) +
    geom_point(alpha=.5, size = 3, aes(colour=as.factor(celltype))) +
    theme(legend.position="none") +
    ggtitle(title)
  return(p) 
  
}

countdata = read.table(inputfile, sep="\t", comment.char = "", header=T, check.names=FALSE, stringsAsFactors=F)
rownames(countdata) = countdata[,1]
countdata = countdata[,-c(1,2)]
colnames(countdata) = str_replace_all(colnames(countdata), ".isoforms.results", "")

counts_dettsnum = data.frame(colSums(countdata >0))
colnames(counts_dettsnum) = c("NumOfTs")
counts_dettsnum$samplename = rownames(counts_dettsnum)
counts_dettsnum = counts_dettsnum[order(counts_dettsnum$samplename),]

  ylab_name = "Number of RSEM detected transcripts"
  title = "Number of RSEM detected transcripts (TPM>0)"
  file_name = "barplot_num_of_ts_rsem"

### draw barplot
g = ggplot(counts_dettsnum, aes(x=samplename,y=NumOfTs)) +
    geom_bar(alpha=0.7, stat="identity") +
    xlab("Sample") + ylab(ylab_name) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle(title)

ggsave(file = paste0(file_name, ".pdf"), plot=g, dpi=100, width=12, height=5)

# Write correlation values to file
outputdata = data.frame(name = counts_dettsnum$samplename, NumOfTs=counts_dettsnum$NumOfTs)
rownames(outputdata) = outputdata$name
write.csv(outputdata, paste0(file_name, ".csv"), quote=FALSE, append=TRUE)


# Printing sessioninfo to standard out
print("Draw plot using TPM counts info:")
sessionInfo()







