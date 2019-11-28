library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

data = read.table("GSE41446.txt",sep="\t",header=TRUE)
datExpr0 = as.data.frame(t(data[, -c(1)]));
names(datExpr0) = data[,1];
rownames(datExpr0) = names(data)[-c(1)];

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)


traitData <- data.frame('SampleId' = names(data)[2:45],"sampleNumber"=1:44,"days" = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7),"surface" = c(1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2))


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitData[1:44,2:4,drop=FALSE], signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
sizeGrWindow(12,9)
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(traitData[1:44,2:4,drop=FALSE]),
main = "Sample dendrogram and trait heatmap")

datExpr = datExpr0
datTraits = traitData

save(datExpr, datTraits, file = "Genes-01-dataInput.RData")
