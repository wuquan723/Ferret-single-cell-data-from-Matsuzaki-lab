library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)
# cluster anslysis---------------------------------------------------------
##This scirpt is for production of Figure 2
##Download and release file into your own file

## Load the datasets by change the PATH to "Ferret_Chromium_E25_P10)
setwd("~/Ferret_sc_data/2020")


E25.data<- Read10X(data.dir = "2020_CR_P367_05_1_E25/outs/filtered_gene_bc_matrices/MusPutFur2.60")
E34_AG.data<-Read10X(data.dir = "2020_CR_P367_01_1_E34GFPP/outs/filtered_gene_bc_matrices/MusPutFur2.60")
E34_neg.data<-Read10X(data.dir = "2020_CR_P367_01_1_E34nega/outs/filtered_gene_bc_matrices/MusPutFur2.60")
E40_neg.data<- Read10X(data.dir = "2020_CR_P318_02_1_E40Neg/outs/filtered_gene_bc_matrices/MusPutFur2.60")
E40_AG.data<- Read10X(data.dir = "2020_CR_P340_01_1_E40GFPP/outs/filtered_gene_bc_matrices/MusPutFur2.60")
P1_AG.data<- Read10X(data.dir = "2020_CR_P400_01_P1GFPpositiveR/outs/filtered_gene_bc_matrices/MusPutFur2.60")
P1_neg.data<- Read10X(data.dir = "2020_CR_P400_02_P1GFPnegR/outs/filtered_gene_bc_matrices/MusPutFur2.60")
P5_AG.data<- Read10X(data.dir = "2020_CR_P382_07_1_P5GFPpositive/outs/filtered_gene_bc_matrices/MusPutFur2.60")
P5_neg.data<- Read10X(data.dir = "2020_CR_P382_06_1_P5GFPnega/outs/filtered_gene_bc_matrices/MusPutFur2.60")
P10_AG.data<- Read10X(data.dir = "2020_CR_P382_08_1_P10GFPpositive/outs/filtered_gene_bc_matrices/MusPutFur2.60")
P10_neg.data<- Read10X(data.dir = "2020_CR_P400_10_P10GFPneg/outs/filtered_gene_bc_matrices/MusPutFur2.60")


E25 <- CreateSeuratObject(E25.data, min.cells = 3, min.features = 200, 
                          project = "10XE25")
E34AG <- CreateSeuratObject(E34_AG.data, min.cells = 3, min.features = 200, 
                            project = "10XE34AG")
E34T <- CreateSeuratObject(E34_neg.data, min.cells = 3, min.features = 200, 
                           project = "10XE34T")
E40AG <- CreateSeuratObject(E40_AG.data, min.cells = 3, min.features = 200, 
                            project = "10XE40AG")
E40T <- CreateSeuratObject(E40_neg.data, min.cells = 3, min.features = 200, 
                           project = "10XE40T")
P1AG <- CreateSeuratObject(P1_AG.data, min.cells = 3, min.features = 200, 
                           project = "10XP1AG")
P1T <- CreateSeuratObject(P1_neg.data, min.cells = 3, min.features = 200, 
                          project = "10XP1T")
P5AG <- CreateSeuratObject(P5_AG.data, min.cells = 3, min.features = 200, 
                           project = "10XP5AG")
P5T <- CreateSeuratObject(P5_neg.data, min.cells = 3, min.features = 200, 
                          project = "10XP5T")
P10AG <- CreateSeuratObject(P10_AG.data, min.cells = 3, min.features = 200, 
                            project = "10XP10AG")
P10T <- CreateSeuratObject(P10_neg.data, min.cells = 3, min.features = 200, 
                           project = "10XP10T")

##option: remove files that are no longer required
rm(list=setdiff(ls(),c("E25","E34AG","E34T","E40AG","E40T","P1AG","P1T","P5AG","P5T","P10AG","P10T")))

Merge <- merge(E25, y = c(E34AG,E34T,E40AG,E40T,P1AG,P1T,P5AG,P5T,P10AG,P10T), 
                add.cell.ids = c("E25","E34AG","E34T","E40AG","E40T","P1AG","P1T","P5AG","P5T","P10AG","P10T"), 
                project = "ferret")

## check cell number
table(Merge$orig.ident)

#10XE25 10XE34AG  10XE34T 10XE40AG  10XE40T 
#3493     3223     2260     1103     1581

#10XP10AG  10XP10T  10XP1AG   10XP1T  10XP5AG   10XP5T 
#3028     3907     2681     3644     2431     2928 

##option: option: remove files that are no longer required
rm(list=setdiff(ls(),"Merge"))

# calculate the percent.mito values.
Merge[["percent.mt"]] <- PercentageFeatureSet(Merge, pattern = "^MT-")
percent.mito <- grep(pattern = "^MT-", x = rownames(GetAssayData(Merge)), value = TRUE) #8 genes

VlnPlot(Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by="orig.ident")

plot1 <- FeatureScatter(Merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))


# We filter out cells that have unique gene counts over 2,500 or less than 200 
# for determining threshold of nGene on GenePlot: always keep within the line nGenevsnUMI
#we don't know Merge mito genes so better to include all cells even with high percent.mito for Merge samples
Merge <- subset(Merge, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)

Merge <- NormalizeData(Merge, normalization.method = "LogNormalize", scale.factor = 10000)

Merge <- FindVariableFeatures(Merge, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

top10 <- head(VariableFeatures(Merge),10)

plot1 <- VariableFeaturePlot(Merge)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1,plot2))

all.genes <- rownames(Merge)
Merge <- ScaleData(Merge, features = all.genes)
Merge <- RunPCA(Merge, features = VariableFeatures(object = Merge))
# Examine and visualize PCA results a few different ways
print(Merge[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Merge, dims = 1:5, reduction = "pca")
DimPlot(Merge, reduction = "pca", dims = c(1, 2))
DimHeatmap(Merge, dims = 1:10, cells = 500, balanced = TRUE)
DimHeatmap(Merge, dims = 11:20, cells = 500, balanced = TRUE)
DimHeatmap(Merge, dims = 21:30, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
Merge <- JackStraw(Merge, num.replicate = 100)
Merge <- ScoreJackStraw(Merge, dims = 1:20)
JackStrawPlot(Merge, dims = 1:20)

ElbowPlot(Merge)

Merge <- FindNeighbors(Merge, dims = 1:20)

Merge <- FindClusters(Merge, resolution = 0.8)
Merge <- RunUMAP(Merge, dims = 1:20)
DimPlot(Merge, reduction = "umap",group.by ="orig.ident")
DimPlot(Merge, reduction = "umap",group.by ="ident", label = TRUE)

table(Idents(Merge))

#clusters
#early RG: 10, 19, 21 (E25)
#mid RG: 7 (E34)
#late RG: 1, 12 ,14
#tRG: 20
#eomes+ IP:  4, 8, 9,15
#Deep Layer neurons: 13 (E25)
#Upper layer neurons:  0, 3, 6, 17
#interneurons: 2
#EP: 23
#OPC: 11
#microglia: 5, 16, 24
#endothelial: 18, 25 (FN1, KDR+ etc), 
#mural cells:  22 (FN1, ACTA2, RGS5, NOTCH3+ etc) (many are common with endothelial genes)
#unknown: 26
new.cluster.ids <- c("F_Upper Layer1","F_lateRG1","F_ITN","F_Upper Layer2","F_IPC1","F_Microglia1",
                     "F_Upper Layer3","F_midRG","F_IPC2","F_IPC3","F_earlyRG1","F_OPC","F_lateRG2",
                     "F_Deep Layer","F_lateRG3","F_Upper Layer4", "F_Microglia2","F_Upper Layer5",
                     "F_Endothelial1","F_earlyRG2","F_tRG","F_earlyRG3","F_mural","F_Ependymal","F_Microglia3","F_Endothelial2","F_unknown")

#assigning cell types identities to clusters before CCA with human data 
#use same cluster name as human one with an additional F at the beginning

names(new.cluster.ids) <- levels(Merge)
Merge <- RenameIdents(Merge, new.cluster.ids)
##Figure 2B and C
DimPlot(Merge, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


##Plot for marker genes in Figure 2E
FeaturePlot(Merge, features = c("PAX6", "EOMES", "TBR1",
                                 "SATB2", "OLIG2", "FOXJ1", "AIF1","GAD1"))
saveRDS(Merge, file = "~/Ferret_sc_data/2020/Merge.rds")
DimPlot(Ferret, reduction = "umap",group.by ="ident", label = TRUE)

# find markers for every cluster compared to all remaining cells, report only the positive ones
Ferretall.markers <- FindAllMarkers(Merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(x = Ferretall.markers,file="Ferretall.markers.txt",sep="\t",row.names=T)
top10 <- Ferretall.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Merge, features = top10$gene) + NoLegend()
