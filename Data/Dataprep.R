## Open libaries and load data ##
detach("package:pd.hugene.1.1.st.v1", unload=TRUE)
detach("package:oligo", unload=TRUE)  # oligo package conflicts with WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# open working directory and load data
setwd("/Users/20172805/Documents/BMT3/OGO comp/Studio") 
data = read.table("Processed data/95varfiltereddata.txt",sep="\t",header=TRUE,row.names=1)

# turn transposed datamatrix into dataframe and transfer row and column names
datExpr0 = as.data.frame(t(data));
names(datExpr0) = rownames(data);
rownames(datExpr0) = names(data);


## Finding outliers ##
# clustering
sampleTree = hclust(dist(datExpr0), method = "average");

# Plot the sample tree
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)


## removing outliers ##
# if necessary #


## creating dataframe with traits data ##
# gender 0 = female, surface 1 = OS, surface 2 = TIO
traitData <- data.frame("gender" = c(0,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0) ,"days" = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7),"surface" = c(1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2))

# creating easier to read sample names
donors = array()
for (i in 1:nrow(traitData)){
  donors[i] = paste('Donor', i,sep = "")
}
rownames(traitData) = donors

# same for datExpr0
samplenames = array()
for (i in 1:nrow(datExpr0)){
  samplenames[i] = paste(rownames(traitData[i,]),"_g",traitData$gender[i],"_d",traitData$days[i],"_s",traitData$surface[i],sep="")
}
rownames(datExpr0) = samplenames


##optional## filtering data based on traits ##
#datExpr3_1 = datExpr[traitData$days == 3 & traitData$surface==1,]
#datExpr3_2 = datExpr[traitData$days == 3 & traitData$surface==2,]
#datExpr7_1 = datExpr[traitData$days == 7 & traitData$surface==1,]
#datExpr7_2 = datExpr[traitData$days == 7 & traitData$surface==2,]

#datExpr3 = datExpr[traitData$days == 3,]
#datExpr7 = datExpr[traitData$days == 7,]

datExprn = datExpr0


## visualising final dataset ##
# Re-cluster samples
sampleTree2 = hclust(dist(datExprn), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitData, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(traitData),
main = "Sample dendrogram and trait heatmap")


## save data to Rdata-file
datExpr = datExprn
datTraits = traitData

save(datExpr, datTraits, file = "Processed data/Genes-dataInput-var95-all.RData")
