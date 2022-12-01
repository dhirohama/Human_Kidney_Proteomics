setwd("~/single_cell/")
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggpubr)

rm(list = ls())
memory.limit(size=56000)  ## set memory


##############################################################################
# data import and data process
##############################################################################

data <- LoadH5Seurat("521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat", assays = "data")

## subset of DKD and LivingDonor(LD)
data <- subset(data, subset = (sampletype == "DKD" | sampletype == "LD"))
data$sampletype <- factor(data$sampletype, levels = c("LD","DKD"))
data <- SetIdent(data, value = data@meta.data$subclass.l1)         ## analysis using subclass.l1



##############################################################################
# UMAP, feature plot, violin plot, bubble plot
##############################################################################

# UMAP
pdf("Umap_KPMP(annotation_DKD_LD).pdf",height=5,width=7)        # with legend
DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE, order = c("Immune","Interstitial","IC","PC","CNT","DCT","TAL","ATL/TAL","DTL","PT","POD","PEC","EC"))			
dev.off()

pdf("Umap_KPMP(annotation_DKD_LD)_1.pdf",height=5,width=5.2)    # no legend
DimPlot(data, reduction = "umap", label = FALSE, repel = TRUE) + NoLegend()			
dev.off()

## Feature plot of MMP7
pdf("feature_MMP7(split.by).pdf",height=5,width=9)
FeaturePlot(data, features = c("MMP7"), split.by = "sampletype") 
dev.off()

## Violin plot of MMP7
pdf("Violin_MMP7(split.by).pdf",height=5,width=15)
VlnPlot(data, features = c("MMP7"), split.by = "sampletype", pt.size = 0)
dev.off()

####  bubble plot
####  In order to Set the standard dot to be the same size for the two figures, then make the figure in the Illustrator
data1 <- subset(data, subset = (sampletype == "LD"))
data2 <- subset(data, subset = (sampletype == "DKD"))

pdf("Bubble_MMP7_LD.pdf",height=5,width=7)
DotPlot(data1, features = c("MMP7"),dot.scale = 6.2, col.min = -0.4, col.max = 2.0) + RotatedAxis() +coord_flip()
dev.off()

pdf("Bubble_MMP7_DKD.pdf",height=5,width=8)
DotPlot(data2, features = c("MMP7"),dot.scale = 7, col.min = -0.4, col.max = 2.0, cols = c("lightgrey","red")) + RotatedAxis() +coord_flip()
dev.off()



##############################################################################
# DEG detection in all 13 clusters
##############################################################################

## PT
PT <- subset(data, idents = "PT")
Idents(PT)=PT@meta.data$sampletype
PT.markers <- FindMarkers(PT, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (PT.markers, file = "PT(subclass.l1).markers.csv")

## CNT
CNT <- subset(data, idents = "CNT")
Idents(CNT)=CNT@meta.data$sampletype
CNT.markers <- FindMarkers(CNT, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (CNT.markers, file = "CNT(subclass.l1).markers.csv")

## PC  
PC <- subset(data, idents = "PC")
Idents(PC)=PC@meta.data$sampletype
PC.markers <- FindMarkers(PC, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (PC.markers, file = "PC(subclass.l1).markers.csv")

## DTL
DTL <- subset(data, idents = "DTL")
Idents(DTL)=DTL@meta.data$sampletype # definition of idents
DTL.markers <- FindMarkers(DTL, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (DTL.markers, file = "DTL(subclass.l1).markers.csv")

## PEC
PEC <- subset(data, idents = "PEC")
Idents(PEC)=PEC@meta.data$sampletype
PEC.markers <- FindMarkers(PEC, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (PEC.markers, file = "PEC(subclass.l1).markers.csv")

## POD
POD <- subset(data, idents = "POD")
Idents(POD)=POD@meta.data$sampletype
POD.markers <- FindMarkers(POD, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (POD.markers, file = "POD(subclass.l1).markers.csv")

## ATL_TAL
ATL_TAL <- subset(data, idents = "ATL/TAL")
Idents(ATL_TAL)=ATL_TAL@meta.data$sampletype
ATL_TAL.markers <- FindMarkers(ATL_TAL, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (ATL_TAL.markers, file = "ATL_TAL(subclass.l1).markers.csv")

## TAL
TAL <- subset(data, idents = "TAL")
Idents(TAL)=TAL@meta.data$sampletype
TAL.markers <- FindMarkers(TAL, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (TAL.markers, file = "TAL(subclass.l1).markers.csv")

## DCT
DCT <- subset(data, idents = "DCT")
Idents(DCT)=DCT@meta.data$sampletype
DCT.markers <- FindMarkers(DCT, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (DCT.markers, file = "DCT(subclass.l1).markers.csv")

## IC
IC <- subset(data, idents = "IC")
Idents(IC)=IC@meta.data$sampletype
IC.markers <- FindMarkers(IC, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (IC.markers, file = "IC(subclass.l1).markers.csv")

## Interstitial
Interstitial <- subset(data, idents = "Interstitial")
Idents(Interstitial)=Interstitial@meta.data$sampletype
Interstitial.markers <- FindMarkers(Interstitial, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (Interstitial.markers, file = "Interstitial(subclass.l1).markers.csv")

## Immune
Immune <- subset(data, idents = "Immune")
Idents(Immune)=Immune@meta.data$sampletype
Immune.markers <- FindMarkers(Immune, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (Immune.markers, file = "Immune(subclass.l1).markers.csv")

## EC
EC <- subset(data, idents = "EC")
Idents(EC)=EC@meta.data$sampletype
EC.markers <- FindMarkers(EC, ident.1 = "DKD", ident.2 = "LD", min.pct = 0.2, logfc.threshold = 0.25)
write.csv (EC.markers, file = "EC(subclass.l1).markers.csv")



##############################################################################
# Box plot of MMP7
##############################################################################

## PT
MMP7_PT <- FetchData(object = PT, vars = c("MMP7"))
MMP7_PT <- MMP7_PT %>%
  mutate(Group = PT@meta.data$sampletype)

compare_means(MMP7~Group, MMP7_PT, method = "wilcox.test", paired = FALSE,
              group.by = NULL)
MMP7_PT$Group <- factor(MMP7_PT$Group, levels = c("LD", "DKD"))

pdf("Box_MMP7_KPMP(PT).pdf",height=5,width=4.5)
p <- ggboxplot(MMP7_PT, x = "Group", y = "MMP7",
               color = "Group",
               add = "jitter",
               add.params = list(size = 0.2)) +
  labs(title="MMP7_KPMP_in_PT",x="Group", y = "Normalized_MMP7")
#  Add p-value
p + stat_compare_means(label.x = 1.3, label.y = 6)
dev.off()

## CNT
MMP7_CNT <- FetchData(object = CNT, vars = c("MMP7"))
MMP7_CNT <- MMP7_CNT %>%
  mutate(Group = CNT@meta.data$sampletype)

compare_means(MMP7~Group, MMP7_CNT, method = "wilcox.test", paired = FALSE,
              group.by = NULL)
MMP7_CNT$Group <- factor(MMP7_CNT$Group, levels = c("LD", "DKD"))

pdf("Box_MMP7_KPMP(CNT).pdf",height=5,width=4.5)
p <- ggboxplot(MMP7_CNT, x = "Group", y = "MMP7",
               color = "Group",
               add = "jitter",
               add.params = list(size = 0.3)) +
  labs(title="MMP7_KPMP_in_CNT",x="Group", y = "Normalized_MMP7")
#  Add p-value
p + stat_compare_means(label.x = 1.3, label.y = 6)
dev.off()

## PC
MMP7_PC <- FetchData(object = PC, vars = c("MMP7"))
MMP7_PC <- MMP7_PC %>%
  mutate(Group = PC@meta.data$sampletype)

compare_means(MMP7~Group, MMP7_PC, method = "wilcox.test", paired = FALSE,
              group.by = NULL)
MMP7_PC$Group <- factor(MMP7_PC$Group, levels = c("LD", "DKD"))

pdf("Box_MMP7_KPMP(PC).pdf",height=5,width=4.5)
p <- ggboxplot(MMP7_PC, x = "Group", y = "MMP7",
               color = "Group",
               add = "jitter",
               add.params = list(size = 0.3)) +
  labs(title="MMP7_KPMP_in_PC",x="Group", y = "Normalized_MMP7")
#  Add p-value
p + stat_compare_means(label.x = 1.3, label.y = 6)
dev.off()
