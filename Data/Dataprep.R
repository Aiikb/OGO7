## Open libaries and load data ##
detach("package:pd.hugene.1.1.st.v1", unload=TRUE)
detach("package:oligo", unload=TRUE)  # oligo package conflicts with WGCNA
library(genefilter)
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()


# open working directory and load data
setwd("/Users/20172805/Documents/BMT3/OGO comp/Studio") 
data = read.table("Processed data/unfiltereddata.txt",sep="\t",header=TRUE,row.names=1)

# turn transposed datamatrix into dataframe and transfer row and column names
datExpr0 = as.data.frame(t(data));
names(datExpr0) = rownames(data);
rownames(datExpr0) = names(data);


## link transcript ids to genenames ##
annot <- read.csv("Raw data/HuGene-1_1-st-v1.na36.hg19.probeset.csv", header=TRUE,
                  skip = 22)
# splitting the geneassignment strings to get gene names
geneAssignment = strsplit(annot$gene_assignment,"//")

# removing probes that have no exon id's
geneNames0 = array()
iOut = array()
for (i in 1:length(geneAssignment)){
  assignment = geneAssignment[[i]]
  if (annot$exon_id[i] == 0){
    iOut = c(iOut,i)
  }
  geneNames0[i] = trimws(assignment[2])
}
iOut = iOut[-1]
geneNames = geneNames0[-iOut]

probes=names(datExpr0)
probes2annot = match(probes, annot$transcript_cluster_id[-iOut])
# check if all probes are annotated
sum(is.na(probes2annot))
# removing unannotated probes
p2aNA <- is.na(probes2annot)
noNAp2a = probes2annot[!p2aNA]
datExprf = datExpr0[!p2aNA]
#probes = names(datExprf)
names(datExprf) = geneNames[noNAp2a]

datExpr0 = datExprf

## creating dataframe with traits data ##
# gender 0 = female, surface 0 = OS, surface 1 = TIO
traitData <- data.frame("gender" = c(0,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0) ,"days" = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7),"surface" = c(0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1))
# creating easier to read sample names
donors = array()
for (i in 1:11){
  donors[i] = paste('Donor', i,sep = "")
}

# add them as rownames
samplenames = array()
for (i in 1:nrow(datExpr0)){
  samplenames[i] = paste(c(donors,donors,donors,donors)[i],"_g",traitData$gender[i],"_d",traitData$days[i],"_s",traitData$surface[i],sep="")
}
rownames(datExpr0) = samplenames
rownames(traitData) = samplenames

##optional## filtering data based on traits ##
#datExpr3_1 = datExpr0[traitData$days == 3 & traitData$surface==1,]
#datExpr3_2 = datExpr0[traitData$days == 3 & traitData$surface==2,]
#datExpr7_1 = datExpr0[traitData$days == 7 & traitData$surface==1,]
#datExpr7_2 = datExpr0[traitData$days == 7 & traitData$surface==2,]

datExpr3 = datExpr0[traitData$days == 3 & traitData$gender == 0,]; traitData3 = traitData[traitData$days == 3 & traitData$gender == 0,]
datExpr7 = datExpr0[traitData$days == 7 & traitData$gender == 0,]; traitData7 = traitData[traitData$days == 7 & traitData$gender == 0,]


datExprn = datExpr3
traitDatan = traitData3

## Finding outliers ##
# clustering
sampleTree = hclust(dist(datExprn));

# Plot the sample tree
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


## removing outliers ##
# if necessary #
#datExprn = datExprn[-c(12,1),]
#traitDatan = traitDatan[-c(12,1),]


data = t(as.matrix(datExprn))

dataF = varFilter(data,var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)

datExprn = as.data.frame(t(dataF))

Genestop100 = names(datExprn)

write(Genestop100,file = "top100day.txt")
write.table(datExprn,file= "top100 day exp.txt")


## visualising final dataset ##
# Re-cluster samples
sampleTree2 = hclust(dist(datExprn))
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitDatan, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(traitDatan),
main = "Sample dendrogram and trait heatmap")



## save data to Rdata-file
datExpr = datExprn
datTraits = traitDatan



save(datExpr, datTraits, file = "Processed data/Genes-dataInput-var50-day3.RData")
