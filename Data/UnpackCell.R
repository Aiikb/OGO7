library(oligo)
library(genefilter)

setwd("/Users/20172805/Documents/BMT3/OGO comp/Studio")               # set working directory to datafolder
untar("Raw data/GSE41446_RAW.tar", exdir = "Raw data/GSE41446/CEL")   # open 
celFiles <- list.files("Raw data/GSE41446/CEL",full=TRUE)     				# open celfiles to celFiles
affyRaw <- read.celfiles(celFiles)				                            # unpack celfiles

eset <- rma(affyRaw)						                                      # normalise data
dataF = varFilter(eset,var.func=IQR, var.cutoff=0.95, filterByQuantile=TRUE)  #filter data by percentile variance

write.exprs(dataF,file="Processed data/95varfiltereddata.txt")				# write to file and log2 transform


