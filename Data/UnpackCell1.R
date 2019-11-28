library(oligo)
library(genefilter)
setwd("/Users/20172805/Documents/BMT3/OGO comp") # set working directory to datafolder
untar("Raw data/GSE41446_RAW.tar", exdir = "GSE41446/CEL")
celFiles <- list.files("GSE41446/CEL",full=TRUE)     				# open celfiles to celFiles
affyRaw <- read.celfiles(celFiles)				# unpack celfiles
eset <- rma(affyRaw)						# normalise data
dataF = varFilter(eset,var.func=IQR, var.cutoff=0.99, filterByQuantile=TRUE)

write.exprs(dataF,file="GSE41446.txt")				# write to file and log2 transform