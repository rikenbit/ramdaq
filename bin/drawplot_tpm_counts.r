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
if (!require("Rtsne")){
    install.packages("Rtsne", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("Rtsne")
}
if (!require("umap")){
    install.packages("umap", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("umap")
}

inputfile <- args[1]
tool_name <- args[2]

### function ###

SDset <- function(x, y){
  data.frame(x=x, y=y, upper = mean(y)+3*sd(y),lower = mean(y)-3*sd(y))
}

create_pca_tsne_umap_mode <- function(countdata, perplexity, mode, local_connectivity=1.0, n_neighbors=15){
  
  tryCatch({
    if (mode=="pca"){
      ### pca
      data_pca = prcomp(t(countdata), scale = FALSE, center = TRUE)
      data_pca_df = as.data.frame(data_pca$x)
      return(list(data=data_pca, data_df=data_pca_df))
    } else if (mode=="tsne") {
      ### tSNE
      data_tsne = Rtsne(t(countdata), dims=2, perplexity=perplexity)
      data_tsne_df = as.data.frame(data_tsne$Y)
      return(list(data=data_tsne, data_df=data_tsne_df))
    } else {
      ### UMAP
      data_umap = umap(t(countdata), local_connectivity=local_connectivity, n_neighbors=n_neighbors)
      data_umap_df = as.data.frame(data_umap$layout)
      return(list(data=data_umap, data_df=data_umap_df))
    }
  }, error = function(e) {
      message("No ", mode, " plot are shown when the number of samples is too small.") 
      print(e)
      return(NULL)
  },finally={})
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

if (tool_name=="rsem"){
  rownames(countdata) = countdata[,1]
  if (ncol(countdata) > 3){
    countdata = countdata[,-c(1,2)]
  } else {
    tmp_rowname = rownames(countdata)
    tmp_colname = colnames(countdata)
    countdata = data.frame(tmp = countdata[,-c(1,2)])
    colnames(countdata) = tmp_colname[3]
    rownames(countdata) = tmp_rowname
  }
  colnames(countdata) = str_replace_all(colnames(countdata), ".genes.results", "")
}

counts_detgenenum = data.frame(colSums(countdata >0))
colnames(counts_detgenenum) = c("NumOfGenes")
counts_detgenenum$samplename = rownames(counts_detgenenum)
counts_detgenenum = counts_detgenenum[order(counts_detgenenum$samplename),]

if (tool_name=="fcounts") {
  ylab_name = "Number of detected genes"
  title = "Number of detected genes (TPM>0)"
  file_name = "barplot_num_of_detectedgene"
} else {
  ylab_name = "Number of RSEM detected genes"
  title = "Number of RSEM detected genes (TPM>0)"
  file_name = "barplot_num_of_gene_rsem"
}

### draw barplot
g = ggplot(counts_detgenenum, aes(x=samplename,y=NumOfGenes)) +
    geom_bar(alpha=0.7, stat="identity") +
    xlab("Sample") + ylab(ylab_name) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle(title)

ggsave(file = paste0(file_name, ".pdf"), plot=g, dpi=100, width=12, height=5)

# Write correlation values to file
outputdata = data.frame(name = counts_detgenenum$samplename, NumOfGenes=counts_detgenenum$NumOfGenes)
rownames(outputdata) = outputdata$name
write.csv(outputdata, paste0(file_name, ".csv"), quote=FALSE, append=TRUE)

### draw pca, tsne, umapplot

if (tool_name=="fcounts"){
  if (ncol(countdata)==1) {
    print("No plot are shown when the number of samples is one.\n") 
  } else {
    countdata_log = log10(countdata+1)
    countdata_pca_data = create_pca_tsne_umap_mode(countdata_log, mode="pca", perplexity=10, local_connectivity=1.0, n_neighbors=15)

    if (!is.null(countdata_pca_data)){
      g1 = ggplot_2D(countdata_pca_data[["data_df"]], "PC1", "PC2", "DRplot_pca_allsample", "")
      ggsave(file = "DRplot_pca_allsample.pdf", plot=g1, dpi=100, width=8, height=7)
      outputdata = countdata_pca_data[["data_df"]][,1:2]
      rownames(outputdata) = colnames(countdata_log)
      write.csv(outputdata, 'DRplot_pca_allsample.csv', quote=FALSE, append=TRUE)
    } 

    ### draw tsneplot
    countdata_tsne_data = create_pca_tsne_umap_mode(countdata_log, mode="tsne", perplexity=10, local_connectivity=1.0, n_neighbors=15)  

    if (!is.null(countdata_tsne_data)){
      g2 = ggplot_2D(countdata_tsne_data[["data_df"]], "V1", "V2", "DRplot_tsne_allsample", "")
      ggsave(file = "DRplot_tsne_allsample.pdf", plot=g2, dpi=100, width=8, height=7)  

      outputdata = countdata_tsne_data[["data_df"]][,1:2]
      rownames(outputdata) = colnames(countdata_log)
      write.csv(outputdata, 'DRplot_tsne_allsample.csv', quote=FALSE, append=TRUE)
    } 

    ### draw umapplot
    countdata_umap_data = create_pca_tsne_umap_mode(countdata_log, mode="umap", perplexity=10, local_connectivity=1.0, n_neighbors=15)  

    if (!is.null(countdata_umap_data)){
      g3 = ggplot_2D(countdata_umap_data[["data_df"]], "V1", "V2", "DRplot_umap_allsample", "")
      ggsave(file = "DRplot_umap_allsample.pdf", plot=g3, dpi=100, width=8, height=7)
      outputdata = countdata_umap_data[["data_df"]][,1:2]
      rownames(outputdata) = colnames(countdata_log)
      write.csv(outputdata, 'DRplot_umap_allsample.csv', quote=FALSE, append=TRUE)
    } 
  }
}


# Printing sessioninfo to standard out
print("Draw plot using TPM counts info:")
sessionInfo()







