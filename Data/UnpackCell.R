library(oligo)
setwd("/Users/20172805/Documents/BMT3/OGO comp/data") # set working directory to datafolder
celFiles <- list.celfiles()     				# open celfiles to celFiles
affyRaw <- read.celfiles(celFiles)				# unpack celfiles
eset <- rma(affyRaw)						# normalise data
write.exprs(eset,file="data.txt")				# write to file and log2 transform


library(lumi)
lumiR("Raw_data_matrix_file.txt", sep = NULL) 