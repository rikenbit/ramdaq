library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(ggrepel)
library(directlabels)
library(entropy)

###########################################################################
### Definition of each variable
###########################################################################

metadata_path = "./Sample_data_sheet.txt"
merged_counts_path = "./merged_SIRVome_data.txt"

#Specify the folder where the plot will be output
outdir = "./output_dir_name"

###########################################################################
### function (Need to run it in advance)
###########################################################################

culc_tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

cv <- function(x) 100*( sd(x)/mean(x))

trim_count_SIRV <- function(data, exclude_id){
  
  data = subset(data, !grepl("ERCC-", gene_name))
  rownames(data) = data$Geneid
  samplename = str_replace_all(colnames(data), "\\.sirv", "")
  colnames(data) = c(samplename)
  
  data = data[,!colnames(data) %in% exclude_id]
  data = data[,colnames(data) != "gene_name"]
  return(data)
}

create_SirvLengthTpm <- function(plotdata, plotpoint=F, outdir, title, output=T){
  
  if(plotpoint){
    data_ends <- plotdata %>% 
      group_by(Geneid) %>% 
      top_n(1, TPM) 
    g = ggplot(plotdata, aes(Length, TPM, colour = GD_RT, group = GD_RT, label = GD_RT)) +
      geom_point() +
      stat_smooth(geom="line", se=FALSE, alpha=0.4, aes(color=GD_RT)) + 
      geom_text_repel(aes(label = Geneid), data = data_ends, size = 3) + 
      theme_linedraw()
    
    if (output){
      ggsave(file = paste0(outdir,"/SirvLengthTpm_w_point_",title,".png"), plot=g, dpi=100, width=12, height=7)
    } else {
      return(g)
    }
    
  } else {
    g = ggplot(plotdata, aes(Length, TPM, colour = GD_RT, group = GD_RT, label = GD_RT)) +
      stat_smooth(geom="line", se=FALSE, alpha=0.8, aes(color=GD_RT)) + 
      theme_linedraw()
    
    if (output){
      ggsave(file = paste0(outdir,"/SirvLengthTpm_wo_point_",title,".png"), plot=g, dpi=100, width=12, height=7)
    } else {
      return(g)
    }
  }
}

###########################################################################
### load metadata & raw counts
###########################################################################

SIRV_metadata = read.table(metadata_path, sep="\t", comment.char = "", header=T, stringsAsFactors=F)
SIRV_counts= read.table(merged_counts_path, sep="\t", comment.char = "", header=T, stringsAsFactors=F)

###########################################################################
### create TPM matrix
###########################################################################

SIRV_counts = trim_count_SIRV(SIRV_counts, "")

SIRV_counts_TPM = culc_tpm(SIRV_counts[,!colnames(SIRV_counts) %in% c("Geneid","Length")], as.numeric(SIRV_counts$Length))

colnames(SIRV_counts_TPM) = paste0(str_split_fixed(colnames(SIRV_counts_TPM),"_", 10)[,1], "_",
                                   str_split_fixed(colnames(SIRV_counts_TPM),"_", 10)[,2], "_",
                                   str_split_fixed(colnames(SIRV_counts_TPM),"_", 10)[,3], "_",
                                   str_split_fixed(colnames(SIRV_counts_TPM),"_", 10)[,4], "_",
                                   str_split_fixed(colnames(SIRV_counts_TPM),"_", 10)[,5], "_",
                                   str_split_fixed(colnames(SIRV_counts_TPM),"_", 10)[,6], "_",
                                   str_split_fixed(colnames(SIRV_counts_TPM),"_", 10)[,7], "_",
                                   str_split_fixed(colnames(SIRV_counts_TPM),"_", 10)[,8], "_",
                                   str_split_fixed(colnames(SIRV_counts_TPM),"_", 10)[,9])
SIRV_counts_TPM = as.data.frame(SIRV_counts_TPM)
SIRV_counts_TPM$Geneid = rownames(SIRV_counts_TPM)

SIRV_TPM_merge = dplyr::left_join(SIRV_counts_TPM, SIRV_counts[,1:2], by = "Geneid")
dim(SIRV_TPM_merge)#84 392

SIRV_TPM_merge_short = SIRV_TPM_merge[1:(84-15),]
SIRV_TPM_merge_long = SIRV_TPM_merge[(84-14):84,]

SIRV_TPM_melt = melt(SIRV_TPM_merge, id=c("Length", "Geneid"))
SIRV_TPM_short_melt = melt(SIRV_TPM_merge_short, id=c("Length", "Geneid"))
SIRV_TPM_long_melt = melt(SIRV_TPM_merge_long, id=c("Length", "Geneid"))

colnames(SIRV_TPM_melt) = c("Length", "Geneid", "Sample", "TPM")
colnames(SIRV_TPM_short_melt) = c("Length", "Geneid", "Sample", "TPM")
colnames(SIRV_TPM_long_melt) = c("Length", "Geneid", "Sample", "TPM")

SIRV_TPM_melt$Group = str_split_fixed(SIRV_TPM_melt$Sample, "_S", 2)[,1]
SIRV_TPM_short_melt$Group = str_split_fixed(SIRV_TPM_short_melt$Sample, "_S", 2)[,1]
SIRV_TPM_long_melt$Group = str_split_fixed(SIRV_TPM_long_melt$Sample, "_S", 2)[,1]

###########################################################################
### create_SirvLengthTpm plot
###########################################################################

SIRV_metadata_trim = SIRV_metadata[,c("Sample_Name", "GD_RT")]

SIRV_TPM_melt = dplyr::left_join(SIRV_TPM_melt, SIRV_metadata_trim, by=c("Sample"="Sample_Name"))
SIRV_TPM_short_melt = dplyr::left_join(SIRV_TPM_short_melt, SIRV_metadata_trim, by=c("Sample"="Sample_Name"))
SIRV_TPM_long_melt = dplyr::left_join(SIRV_TPM_long_melt, SIRV_metadata_trim, by=c("Sample"="Sample_Name"))

create_SirvLengthTpm(SIRV_TPM_melt, F, outdir, "all")
create_SirvLengthTpm(SIRV_TPM_short_melt, F, outdir, "short")
create_SirvLengthTpm(SIRV_TPM_long_melt, F, outdir, "long")

create_SirvLengthTpm(SIRV_TPM_melt, T, outdir, "all")
create_SirvLengthTpm(SIRV_TPM_short_melt, T, outdir, "short")
create_SirvLengthTpm(SIRV_TPM_long_melt, T, outdir, "long")

###########################################################################
### TPM mean & sd
###########################################################################

SIRV_TPM_plotdata = melt(SIRV_TPM_merge, id=c("Length", "Geneid"))
SIRV_TPM_plotdata = SIRV_TPM_plotdata[,c(2,3,4)]
colnames(SIRV_TPM_plotdata) = c( "Geneid", "Sample", "TPM")
SIRV_TPM_plotdata$TPM_log = log10(SIRV_TPM_plotdata$TPM+1)

SIRV_TPM_plotdata = dplyr::left_join(SIRV_TPM_plotdata, SIRV_metadata_trim, by=c("Sample"="Sample_Name"))

SIRV_TPM_summarise = SIRV_TPM_plotdata %>% 
  group_by(Sample) %>%
  summarise_each(funs(mean, sd, cv), -c(Geneid, GD_RT))
SIRV_TPM_summarise = dplyr::left_join(SIRV_TPM_summarise, SIRV_metadata_trim, by=c("Sample"="Sample_Name"))

GD_RT_summarise = SIRV_TPM_plotdata %>% 
  group_by(GD_RT) %>%
  summarise_each(funs(mean, sd, cv), -c(Sample, Geneid))
GD_RT_summarise

data_ends <- SIRV_TPM_plotdata %>% 
  group_by(Geneid) %>% 
  top_n(1, TPM) 

g = ggplot(SIRV_TPM_summarise, aes(TPM_log_mean, TPM_log_sd, colour = GD_RT)) +
  geom_point() +
  geom_point(data=GD_RT_summarise, aes(fill=GD_RT), size=4, shape=24) +
  geom_text_repel(aes(label = GD_RT), data = GD_RT_summarise, size = 3, box.padding = 2, max.overlaps = Inf) + 
  theme_linedraw() 
ggsave(file = paste0(outdir,"/SirvMeanSdLog.png"), plot=g, dpi=100, width=10, height=7)

###########################################################################
### plot entropy
###########################################################################

### unit=log
data_entropy = apply( SIRV_counts_TPM[,-391], 2, entropy.empirical)
data_entropy_set = data.frame(Sample= names(data_entropy), Entropy=data_entropy)
data_entropy_set = dplyr::left_join(data_entropy_set, SIRV_metadata_trim, by=c("Sample"="Sample_Name"))

g = ggplot(data_entropy_set, aes(x=GD_RT, y=Entropy, fill=GD_RT)) + 
  geom_violin()+
  geom_jitter(shape=16) +
  geom_hline(yintercept=4.430817, linetype="dashed", color = "red") + 
  theme_linedraw()
ggsave(file = paste0(outdir,"/SirvTpmEntropy.png"), plot=g, dpi=100, width=12, height=5)












