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

### function

calc_Tcenter <- function(bed.row){
  
  if (bed.row[10]==1){
    new_position = as.integer(bed.row[2] + (bed.row[3]-bed.row[2])/2)
    return(new_position)
  }
  
  t_Length_set = as.integer(unlist(str_split(bed.row[11], ",")))
  t_Length_set = t_Length_set[!is.na(t_Length_set)]
  t_length = as.integer(sum(t_Length_set)/2)
  
  t_Position_set = as.integer(unlist(str_split(bed.row[12], ",")))
  t_Position_set = t_Position_set[!is.na(t_Position_set)]
  
  for(i in 1:length(t_Length_set)){
    S_i = sum(t_Length_set[1:i]) 
    
    if(S_i >= t_length){
      k = i
      s = sum(t_Length_set[1:(i-1)]) 
      new_position = as.integer(bed.row[2] + t_Position_set[k] + (t_length-s))
      return(new_position)
    }
  }
}

inputfile <- args[1]
bed.raw = read.table(inputfile, sep="\t", header=F, stringsAsFactors=F)

bed.adjusted = bed.raw
for(i in 1:nrow(bed.raw)){
  if(bed.raw[i,7] == bed.raw[i,8]){
    bed.adjusted[i,7] = calc_Tcenter(bed.raw[i,])
    bed.adjusted[i,8] = bed.adjusted[i,7]
  }
}

write.table(bed.adjusted, 'adjusted.bed', sep="\t", col.names=F, row.names=F, quote=F)

# Printing sessioninfo to standard out
print("adjusted bed script info:")
sessionInfo()







