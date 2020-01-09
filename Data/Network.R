## Open libaries and load data ##
detach("package:pd.hugene.1.1.st.v1", unload=TRUE)
detach("package:oligo", unload=TRUE)  # oligo package conflicts with WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# open working directory and load data
setwd("/Users/20172805/Documents/BMT3/OGO comp/Studio") 
lnames = load(file = "Processed data/Genes-dataInput-var50-day7.RData")
lnames


## choosing a soft-threshold power ## in case no power below 15 is found see https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
# a set of soft-thresholding powers to try
powers = c(1:20)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
# Scale-free topology fit index as a function of the soft-thresholding power
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


## Creating the network##
net = blockwiseModules(datExpr, maxBlockSize = 10000, power = 9,
TOMType = "unsigned", minModuleSize = 20,
reassignThreshold = 0, mergeCutHeight = 0.15,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "Nets/GeneTOM",
verbose = 3)
table(net$colors)

recutnet = net


#optional# recut modules
minModuleSize = min(20, ncol(datExpr)/2 )
recutnet = recutBlockwiseTrees(
  datExpr,
  net$goodSamples, net$goodGenes,
  net$blocks,
  net$TOMFiles,
  net$dendrograms,
  corType = "pearson",
  networkType = "unsigned",
  deepSplit = 2,
  detectCutHeight = 0.995, minModuleSize = min(20, ncol(datExpr)/2 ),
  maxCoreScatter = NULL, minGap = NULL,
  maxAbsCoreScatter = NULL, minAbsGap = NULL,
  minSplitHeight = NULL, minAbsSplitHeight = NULL,
  
  useBranchEigennodeDissim = FALSE,
  minBranchEigennodeDissim = mergeCutHeight,
  
  pamStage = TRUE, pamRespectsDendro = TRUE,
  minCoreKME = 0.5, minCoreKMESize = minModuleSize/3,
  minKMEtoStay = 0.3,
  reassignThreshold = 1e-6,
  mergeCutHeight = 0.15, impute = TRUE,
  trapErrors = FALSE, numericLabels = FALSE,
  verbose = 0, indent = 0)
table(recutnet$colors)


# Convert labels to colors for plotting
mergedColors = labels2colors(recutnet$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)


## save the network to Rdata-file##
moduleLabels = recutnet$colors
moduleColors = labels2colors(recutnet$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

save(MEs, moduleLabels, moduleColors, geneTree,
file = "Nets/Topography-networkConstruction-auto 50var day7.RData")
