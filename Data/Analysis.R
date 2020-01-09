## Open libaries and load data ##
detach("package:pd.hugene.1.1.st.v1", unload=TRUE)
detach("package:oligo", unload=TRUE)  # oligo package conflicts with WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# open working directory and load data
setwd("/Users/20172805/Documents/BMT3/OGO comp/Studio") # set working directory to datafolder
lnames = load(file = "Processed data/Genes-dataInput-50var-day7.RData")
lnames
lnames = load(file = "Nets/Topography-networkConstruction-auto 50var day7.RData");
lnames

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);


## Correlation module eigengenes with traits##
# Recalculate Module eigengenes with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
# Calculate correlation with traits and p-values
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Create matrix with correlations between modules and traits and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), " (",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor);
par(mar = c(0.5, 9, 2, 2));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = datTraits,
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))


## Correlation genes with traits and modules ##
# Calculate gene-module membership and p values
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# Calculate gene-trait significance and p values
geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

GS_names = array()
pGS_names = array()
for (i in 1:length(datTraits)){
  GS_names[i] = paste("GS.", names(datTraits)[i], sep="")
  pGS_names[i] = paste("p.GS.", names(datTraits)[i], sep="")
}
names(geneTraitSignificance) = GS_names
names(GSPvalue) = pGS_names


# Plot chosen module membership vs. gene significance
#sizeGrWindow(7, 7);
par(mfrow = c(1,1));
module = "tan"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance$GS.surface[moduleGenes]),
corFnc = "cor", corOptions = "use = 'p'",xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for surface",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,abline = TRUE, col = "black")
fit <- lm(abs(geneModuleMembership[moduleGenes, column]) ~ abs(geneTraitSignificance$GS.surface[moduleGenes]))
R2 = summary(fit)$adj.r.squared
legend("topright", bty="n", legend=paste("R2 is",format(R2, digits=4)))



## create summary of geneinfo ##
# creating dataframe with genenames, module colors, genetraitsignificance
geneInfo0 = data.frame(geneName = names(datExpr),
                       moduleColor = moduleColors)


for (trait in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[, trait],
                         GSPvalue[, trait]);
  names(geneInfo0) = c(oldNames, names(geneTraitSignificance)[trait], names(GSPvalue)[trait])
}

# Order modules by their significance for days
modOrder = order(-abs(moduleTraitCor[,3]));
# Add module membership information columns with the most significant first
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.surface));
geneInfo = geneInfo0[geneOrder, ]


##create summary of modinfo ##
# create dataframe with module-trait correlation and p values
modInfo0 = data.frame(moduleTraitCor[modOrder[1],],moduleTraitPvalue[modOrder[1],])
names(modInfo0) = c(paste("MTC.", modNames[modOrder[1]], sep=""),
                    paste("p.MTC.", modNames[modOrder[1]], sep=""))
for (mod in 2:nrow(moduleTraitCor))
{
  oldNames = names(modInfo0)
  modInfo0 = data.frame(modInfo0,moduleTraitCor[modOrder[mod],],moduleTraitPvalue[modOrder[mod],]);
  names(modInfo0) = c(oldNames,paste("MTC.", modNames[modOrder[mod]], sep=""),
                       paste("p.MTC.", modNames[modOrder[mod]], sep=""))
}
# transpose dataframe
modInfo = as.data.frame(t(modInfo0))


## saving data ##
write.csv(geneInfo, file = "Analysis/geneInfo_50var_day7.csv")
write.csv(modInfo, file = "Analysis/modInfo_50var_day7.csv")

