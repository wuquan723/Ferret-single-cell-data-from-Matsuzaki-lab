library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)
library(data.table)
## Load the datasets by change the PATH to "Merge.rds")
setwd("~/Ferret_sc_data/2023_R_script_data_for_share")
Merge <- readRDS(file = "~/Ferret_sc_data/2023_R_script_data_for_share/Merge.rds")

Ferret_NP<- subset(Merge, idents = c("F_lateRG1","F_IPC1","F_midRG","F_IPC2","F_IPC3","F_earlyRG1","F_OPC","F_lateRG2",
                                     "F_lateRG3","F_earlyRG2","F_tRG","F_earlyRG3","F_Ependymal"))


###Go to the website (https://cells.ucsc.edu/?ds=cortex-dev#), you download the exprMatrix.tsv.gz and meta.tsv files.
##from Info&Download
file =gzfile('exprMatrix.tsv.gz','rt')   
expr.data <- read.table(file,header=T,row.names = 1)
expr.data <- fread("exprMatrix.tsv.gz")
meta <- read.csv('meta.tsv',sep="\t")
genes = expr.data[,1][[1]]
genes = gsub(".+[|]", "", genes)
expr.data = data.frame(expr.data[,-1], row.names=genes)

rownames(meta)<-meta$Cell
Human_Sci <- CreateSeuratObject(counts = expr.data, min.cells = 3, min.features = 200, 
                                project = "H_Sci", meta.data = meta)

Human_Sci <- NormalizeData(Human_Sci, normalization.method = "LogNormalize", scale.factor = 10000)

Human_Sci<- FindVariableFeatures(Human_Sci, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
all.genes <- rownames(Human_Sci)

Human_Sci <- ScaleData(Human_Sci, features = all.genes)

Human_Sci <- RunPCA(Human_Sci, features = VariableFeatures(object = Human_Sci))

DimPlot(Human_Sci, reduction = "pca", dims = c(2, 3))
DimHeatmap(Human_Sci, dims = 1:10, cells = 500, balanced = TRUE)
DimHeatmap(Human_Sci, dims = 11:20, cells = 500, balanced = TRUE)

Human_Sci <- FindNeighbors(Human_Sci, dims = 1:20)

Human_Sci <- FindClusters(Human_Sci, resolution = 0.8)
Human_Sci <- RunUMAP(Human_Sci, dims = 1:20)
DimPlot(Human_Sci, reduction = "umap",group.by ="WGCNAcluster")
DimPlot(Human_Sci, reduction = "umap",group.by ="WGCNAcluster") + facet_wrap(~WGCNAcluster)


##Rename work identity
Idents(Human_Sci) <- Human_Sci@meta.data$WGCNAcluster
table(Idents(Human_Sci))
#remove low quality cluster and non-cortical cells (absent in our ferret samples):
Human_Sci_subset <- subset(x = Human_Sci, idents = c("EN-PFC1","EN-PFC2","MGE-IPC3","MGE-RG1","EN-PFC3","EN-V1-1","MGE-IPC1","MGE-IPC2","MGE-RG2",
                                                     "EN-V1-2","EN-V1-3","Glyc","IN-STR","MGE-div","nIN1","nIN2","nIN3","nIN4","nIN5"), invert = TRUE)

DimPlot(Human_Sci_subset, reduction = "umap",group.by ="ident")


#### (2022 January 21) Subset NP Science Human ####
levels(Human_Sci)
Human_Sci_NP <- subset(Human_Sci, idents = c("RG-early","RG-div1","RG-div2","vRG","tRG","oRG","OPC",
                                             "Astrocyte","IPC-div1","IPC-div2","IPC-nEN1","IPC-nEN2","IPC-nEN3"))
#926 cells


table(Human_Sci_NP$WGCNAcluster)
Human_Sci$orig.ident <- Human_Sci$Age_in_Weeks

#### CCA Human dataset from Science paper and Ferret NP ####
anchors_Sci <- FindIntegrationAnchors(object.list =list(Human_Sci,Ferret_NP), dims = 1:30,
                                      k.anchor = 5,
                                      k.filter = 200,
                                      k.score = 30,
                                      max.features = 200,
                                      nn.method = "annoy",
                                      n.trees = 50,
                                      eps = 0)

integrated_Sci <- IntegrateData(anchorset = anchors_Sci, dims = 1:30)
#Merging dataset 1 into 2
#Extracting anchors for merged samples
#Finding integration vectors
#Integrating data

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(integrated_Sci) <- "integrated"

# Run the standard workflow for visualization and clustering
integrated_Sci <- ScaleData(integrated_Sci, verbose = FALSE)
integrated_Sci <- RunPCA(integrated_Sci, npcs = 30, verbose = FALSE)

integrated_Sci <- RunUMAP(integrated_Sci, reduction = "pca", dims = 1:30)
integrated_Sci <- FindNeighbors(integrated_Sci, dims = 1:30)
integrated_Sci <- FindClusters(integrated_Sci, resolution = 0.8)
DimPlot(integrated_Sci, reduction = "umap",group.by ="ident", label = TRUE)+ NoLegend()
DimPlot(integrated_Sci, reduction = "umap", group.by = "ident", label = TRUE)
DimPlot(integrated_Sci, reduction = "umap", group.by = "orig.ident", label = TRUE)
DimPlot(integrated_Sci, reduction = "umap", group.by = "orig.ident") +facet_wrap(~orig.ident, nrow = 5)+ NoLegend()


##save dataset for next step##

saveRDS(Human_Sci_NP, file = "~/Human_Sci_NP.rds")
saveRDS(Ferret_NP, file = "~/Ferret_NP.rds")
