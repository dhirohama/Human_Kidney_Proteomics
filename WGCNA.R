setwd("~/WGCNA/")

# clear R's brain
rm(list = ls())

library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
library(flashClust)
library(gplots)
library(tidyverse)
library(tidyverse, warn.conflicts=F, quietly=T)



##############################################################################
# data import and data process
##############################################################################

# matrix data import
dataExpr <- read.csv("Proteomics_expression.csv", header = TRUE, row.names = 1, fileEncoding="UTF-8-BOM")
## phenotype data import
traitData <- read.csv("Proteomics_traitData.csv", header = TRUE, row.names = 1, fileEncoding="UTF-8-BOM") 

dataExpr <- log2(dataExpr)
dataExpr <- scale(dataExpr)
dataExpr <- data.frame(dataExpr)
dataExpr <- as.data.frame(t(dataExpr))

gsg = goodSamplesGenes(dataExpr, verbose = 3)   ##  Flagging proteins and samples with too many missing values
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
table(rownames(traitData) == rownames(dataExpr))  # TRUE



##############################################################################
# Establishment of sample network, outliers check, and find Co-expression modules
##############################################################################

A = adjacency(t(dataExpr), type = "signed")   
k = as.numeric(apply(A, 2, sum)) - 1
Z.k = scale(k)
thresholdZ.k = -5  # no outliers were found

sampleTree = hclust(dist(dataExpr), method = "average")
traitData <- select(traitData, "eGFR","Groups","DKD","HTN","BMI","Age","Gender","Race","GS","Fibrosis")
traitColors = data.frame(numbers2colors(traitData, signed = FALSE))
dimnames(traitColors)[[2]] = paste(names(traitData), "C", sep = "")

type = "signed" 
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)

pdf("Scale_independence.pdf",width=4, height=5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red",abline(h=0.9,col="red"))
dev.off()

pdf("Mean_connectivity.pdf",width=4, height=5)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     col="red")
dev.off()



##############################################################################
# WGCNA
##############################################################################

softPower = 10   ## power was chosen based on figure above
adjacency = adjacency(dataExpr, power = softPower, type="signed")
TOM = TOMsimilarity(adjacency,TOMType = "signed")
dissTOM = 1-TOM

k=as.vector(apply(adjacency, 2, sum, na.rm=T))  
par(mfrow = c(1,1))
pdf("histK.pdf",width=6, height=5)
hist(k)
dev.off()

pdf("Scalefree.pdf",width=6, height=5)     
scaleFreePlot(k, main="Check scale free topology\n")
dev.off()

geneTree = flashClust(as.dist(dissTOM), method = "average")       # make geneTree

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid", 
                            deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
print(table(dynamicMods))  # confirm the number of modules (=9 modules)

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# merge module based on MEDissThres = 0.25
MEList <- moduleEigengenes(dataExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
MEDissThres = 0.25

# comparison before and after merge
merge <- mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
pdf("sampleDendro_merged.pdf",width=6, height=5)
plotDendroAndColors(geneTree, mergedColors, "Merged dynamic",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()

# sampleTreeHeatmap
plotTOM = dissTOM^7
moduleColors = dynamicColors
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
png("sampleTreeHeatmap.png", width = 640, height = 640)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes", col=myheatcol)
dev.off()

MEs <- mergedMEs
dynamicColors <- mergedColors
table(dynamicColors)

# Plot Eigengenenetwork
# Add eGFR and Fibrosis to existing module eigengenes
eGFR = as.data.frame(traitData$eGFR)
names(eGFR) = "eGFR"
Fibrosis = as.data.frame(traitData$Fibrosis)
names(Fibrosis) = "Fibrosis"
MET = orderMEs(cbind(MEs, eGFR, Fibrosis))

pdf("EigengeneNetworks(heatmap_with_clinical_trait).pdf",width=6, height=5)          # heatmap
plotEigengeneNetworks(MET, "", plotDendrograms = FALSE,marHeatmap = c(6, 
                                                                      6, 1, 2), cex.lab = 0.8, xLabelsAngle = 90) 
dev.off()

pdf("EigengeneNetworks(dendro_with_clinical_trait(eGFR,Fibrosis)).pdf",width=6, height=5)      # dengrogram
plotEigengeneNetworks(MET, "", plotHeatmaps = FALSE,marDendro = c(2, 4, 1, 2), cex.lab = 0.8) 
dev.off()

## further add DKD, BMI, Age, Gender, Race, GS, HTN
DKD = as.data.frame(traitData$DKD)
names(DKD) = "DKD"
BMI = as.data.frame(traitData$BMI)
names(BMI) = "BMI"
Age = as.data.frame(traitData$Age)
names(Age) = "Age"
Gender = as.data.frame(traitData$Gender)
names(Gender) = "Gender"
Race = as.data.frame(traitData$Race)
names(Race) = "Race"
GS = as.data.frame(traitData$GS)
names(GS) = "GS"
HTN = as.data.frame(traitData$HTN)
names(HTN) = "HTN"
MET = orderMEs(cbind(MET, DKD, BMI, Age, Gender, Race, GS, HTN))

pdf("EigengeneNetworks(dendro_with_clinical_trait).pdf",width=6, height=5)      # dengrogram
plotEigengeneNetworks(MET, "", plotHeatmaps = FALSE,marDendro = c(2, 4, 1, 2), cex.lab = 0.8) 
dev.off()


#-----module relationship
#Define number of genes and samples
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
# make correlation between module and trait
moduleTraitCor = cor(MEs, traitData, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)

pdf("modulerelationship_heatmap.pdf",width=12, height=12)
par(mar = c(5, 10, 5, 5))
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(traitData),
               yLabels= names(MEs),　　
               ySymbols= names(MEs),　　
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 1.4,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))
dev.off()



##############################################################################
# pick up module proteins from WGCNA
##############################################################################

#---module gene info
annot = read.csv(file = "SomarmerID_EntrezGeneSymbol.csv", header = TRUE, fileEncoding="UTF-8-BOM");
dim(annot)
names(annot)
probes = names(dataExpr)
probes2annot = match(probes, annot$geneid)  # The following is the number or probes without annotation:
sum(is.na(probes2annot)) # Should return 0.

# Create the starting data frame
geneInfo = data.frame(geneid = probes,
                       geneSymbol = annot$EntrezGeneSymbol[probes2annot],
                       LocusLinkID = annot$EntrezGeneID[probes2annot],
                       moduleColor = dynamicColors)

write.csv(geneInfo, file = "geneInfo_all_module.csv")
