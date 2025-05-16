# Loading the libraries for analysis ####
#Packages for sc-RNA seq - seurat workflow
library(devtools)
library(R.utils)
library(utils)
library(ggplot2)
library(reshape2)
library(readxl)
library(openxlsx)
library(tidyverse)
library(RColorBrewer) 
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(pkgbuild)
library(presto)
library(glmGamPoi)
library(Seurat)
library(scDblFinder)
library(CellChat)
library(hdf5r)
library(monocle3)
library(SeuratWrappers)
# Packages for bulk/pseudobulk RNA seq and enrichment (EnrichR)
library(DESeq2)
library(ggrepel)
library(openxlsx)
library(GEOquery)
library(limma)
library(umap)
library(clusterProfiler)
library(org.Mm.eg.db)
library("ReactomePA")
library(stringr)
library(ggplot2)
library(tidyverse)
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023","Reactome_2022")
library(scales)
library(cowplot)

 future::plan("multisession", workers = 5) # do parallel
 options(future.globals.maxSize = 48000 * 1024^2) #increases the max global export for future expression (if you need to)

# PRE-INTEGRATION: Loading in the single-cell wound healing data set and performing individual data QC ####
KOFD0.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/KOFD0/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)
KOFD3.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/KOFD3/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)
KOFD7.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/KOFD7/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)

KOMD0.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/KOMD0/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)
KOMD3.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/KOMD3/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)
KOMD7.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/KOMD7/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)

WTFD0.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/WTFD0/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)
WTFD3.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/WTFD3/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)
WTFD7.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/WTFD7/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)

WTMD0.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/WTMD0/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)
WTMD3.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/WTMD3/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)
WTMD7.obj<- Read10X_h5("C:/Users/chenh/Desktop/Novogene Sequencing Folders/Novogene sc-RNA seq/Result_X202SC23104597-Z01-F004_mm10/Result_X202SC23104597-Z01-F004_mm10/1.Data_process/WTMD7/filtered_feature_bc_matrix.h5",use.names = TRUE,unique.features = TRUE)

### Creating seuratObject - preQC
KOFD0.SO <- CreateSeuratObject(counts = KOFD0.obj, min.cells = 5, min.features = 200)
KOFD3.SO <- CreateSeuratObject(counts = KOFD3.obj, min.cells = 5, min.features = 200)
KOFD7.SO <- CreateSeuratObject(counts = KOFD7.obj, min.cells = 5, min.features = 200)

KOMD0.SO <- CreateSeuratObject(counts = KOMD0.obj, min.cells = 5, min.features = 200)
KOMD3.SO <- CreateSeuratObject(counts = KOMD3.obj, min.cells = 5, min.features = 200)
KOMD7.SO <- CreateSeuratObject(counts = KOMD7.obj, min.cells = 5, min.features = 200)

WTFD0.SO <- CreateSeuratObject(counts = WTFD0.obj, min.cells = 5, min.features = 200)
WTFD3.SO <- CreateSeuratObject(counts = WTFD3.obj, min.cells = 5, min.features = 200)
WTFD7.SO <- CreateSeuratObject(counts = WTFD7.obj, min.cells = 5, min.features = 200)

WTMD0.SO <- CreateSeuratObject(counts = WTMD0.obj, min.cells = 5, min.features = 200)
WTMD3.SO <- CreateSeuratObject(counts = WTMD3.obj, min.cells = 5, min.features = 200)
WTMD7.SO <- CreateSeuratObject(counts = WTMD7.obj, min.cells = 5, min.features = 200)

# Identifying each seurat object's name
#KO
Idents(KOFD0.SO) <- "KOFD0"
Idents(KOFD3.SO) <- "KOFD3"
Idents(KOFD7.SO) <- "KOFD7"

Idents(KOMD0.SO) <- "KOMD0"
Idents(KOMD3.SO) <- "KOMD3"
Idents(KOMD7.SO) <- "KOMD7"

#WT
Idents(WTFD0.SO) <- "WTFD0"
Idents(WTFD3.SO) <- "WTFD3"
Idents(WTFD7.SO) <- "WTFD7"

Idents(WTMD0.SO) <- "WTMD0"
Idents(WTMD3.SO) <- "WTMD3"
Idents(WTMD7.SO) <- "WTMD7"

# Calculate the percentage of mitochondrial genes in each cell and assign it as a metadata variable
#KO
KOFD0.SO[["percent.mt"]] <- PercentageFeatureSet(KOFD0.SO, pattern = "^mt-")
KOFD3.SO[["percent.mt"]] <- PercentageFeatureSet(KOFD3.SO, pattern = "^mt-")
KOFD7.SO[["percent.mt"]] <- PercentageFeatureSet(KOFD7.SO, pattern = "^mt-")

KOMD0.SO[["percent.mt"]] <- PercentageFeatureSet(KOMD0.SO, pattern = "^mt-")
KOMD3.SO[["percent.mt"]] <- PercentageFeatureSet(KOMD3.SO, pattern = "^mt-")
KOMD7.SO[["percent.mt"]] <- PercentageFeatureSet(KOMD7.SO, pattern = "^mt-")

#WT
WTFD0.SO[["percent.mt"]] <- PercentageFeatureSet(WTFD0.SO, pattern = "^mt-")
WTFD3.SO[["percent.mt"]] <- PercentageFeatureSet(WTFD3.SO, pattern = "^mt-")
WTFD7.SO[["percent.mt"]] <- PercentageFeatureSet(WTFD7.SO, pattern = "^mt-")

WTMD0.SO[["percent.mt"]] <- PercentageFeatureSet(WTMD0.SO, pattern = "^mt-")
WTMD3.SO[["percent.mt"]] <- PercentageFeatureSet(WTMD3.SO, pattern = "^mt-")
WTMD7.SO[["percent.mt"]] <- PercentageFeatureSet(WTMD7.SO, pattern = "^mt-")

#cleaning up the metadata
#Create a new metadata column called "original" with the contents of "orig.ident"
KOFD0.SO$orig.identity <- KOFD0.SO$orig.ident
KOFD3.SO$orig.identity <- KOFD3.SO$orig.ident
KOFD7.SO$orig.identity <- KOFD7.SO$orig.ident

KOMD0.SO$orig.identity <- KOMD0.SO$orig.ident
KOMD3.SO$orig.identity <- KOMD3.SO$orig.ident
KOMD7.SO$orig.identity <- KOMD7.SO$orig.ident

WTFD0.SO$orig.identity <- WTFD0.SO$orig.ident
WTFD3.SO$orig.identity <- WTFD3.SO$orig.ident
WTFD7.SO$orig.identity <- WTFD7.SO$orig.ident

WTMD0.SO$orig.identity <- WTMD0.SO$orig.ident
WTMD3.SO$orig.identity <- WTMD3.SO$orig.ident
WTMD7.SO$orig.identity <- WTMD7.SO$orig.ident

# Remove the "orig.ident" metadata column
KOFD0.SO$orig.identity <- Idents(KOFD0.SO)
KOFD0.SO$orig.ident <- NULL
KOFD3.SO$orig.identity <- Idents(KOFD3.SO)
KOFD3.SO$orig.ident <- NULL
KOFD7.SO$orig.identity <- Idents(KOFD7.SO)
KOFD7.SO$orig.ident <- NULL

KOMD0.SO$orig.identity <- Idents(KOMD0.SO)
KOMD0.SO$orig.ident <- NULL
KOMD3.SO$orig.identity <- Idents(KOMD3.SO)
KOMD3.SO$orig.ident <- NULL
KOMD7.SO$orig.identity <- Idents(KOMD7.SO)
KOMD7.SO$orig.ident <- NULL

WTFD0.SO$orig.identity <- Idents(WTFD0.SO)
WTFD0.SO$orig.ident <- NULL
WTFD3.SO$orig.identity <- Idents(WTFD3.SO)
WTFD3.SO$orig.ident <- NULL
WTFD7.SO$orig.identity <- Idents(WTFD7.SO)
WTFD7.SO$orig.ident <- NULL

WTMD0.SO$orig.identity <- Idents(WTMD0.SO)
WTMD0.SO$orig.ident <- NULL
WTMD3.SO$orig.identity <- Idents(WTMD3.SO)
WTMD3.SO$orig.ident <- NULL
WTMD7.SO$orig.identity <- Idents(WTMD7.SO)
WTMD7.SO$orig.ident <- NULL

#saving the pre-QC Seurat Objects 
#saving the KO
saveRDS(
  object = KOFD0.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/KOFD0.SO.Rds")
saveRDS(
  object = KOFD3.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/KOFD3.SO.Rds")
saveRDS(
  object = KOFD7.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/KOFD7.SO.Rds")

saveRDS(
  object = KOMD0.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/KOMD0.SO.Rds")
saveRDS(
  object = KOMD3.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/KOMD3.SO.Rds")
saveRDS(
  object = KOMD7.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/KOMD7.SO.Rds")

#savng the WT
saveRDS(
  object = WTFD0.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/WTFD0.SO.Rds")
saveRDS(
  object = WTFD3.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/WTFD3.SO.Rds")
saveRDS(
  object = WTFD7.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/WTFD7.SO.Rds")

saveRDS(
  object = WTMD0.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/WTMD0.SO.Rds")
saveRDS(
  object = WTMD3.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/WTMD3.SO.Rds")
saveRDS(
  object = WTMD7.SO,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Pre-QC Seurat Object/WTMD7.SO.Rds")


# Visualize the distribution of detected genes, numbers of RNA and the mitochondrial percentage in all cells (Before moving on to step 2.20...create scatter plots for each of the Seurat Objects)
#KOFD0
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.KOFD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOFD0.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(KOFD0.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(KOFD0.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.KOF0.tiff", units="in", width=30, height=5, res=300)
plot1 + plot2
dev.off()
#KOFD3
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.KOFD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOFD3.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot3 <- FeatureScatter(KOFD3.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(KOFD3.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.KOFD3.tiff", units="in", width=30, height=5, res=300)
plot3 + plot4
dev.off()

#KOFD7
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.KOFD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOFD7.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot5 <- FeatureScatter(KOFD7.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(KOFD7.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.KOFD7.tiff", units="in", width=30, height=5, res=300)
plot5 + plot6
dev.off()

#KOMD0
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.KOMD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOMD0.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(KOMD0.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(KOMD0.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.KOM0.tiff", units="in", width=30, height=5, res=300)
plot1 + plot2
dev.off()

#KOMD3
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.KOMD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOMD3.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot3 <- FeatureScatter(KOMD3.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(KOMD3.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.KOMD3.tiff", units="in", width=30, height=5, res=300)
plot3 + plot4
dev.off()

#KOMD7
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.KOMD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOMD7.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot5 <- FeatureScatter(KOMD7.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(KOMD7.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.KOMD7.tiff", units="in", width=30, height=5, res=300)
plot5 + plot6
dev.off()

#WTFD0
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.WTFD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTFD0.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(WTFD0.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WTFD0.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.WTF0.tiff", units="in", width=30, height=5, res=300)
plot1 + plot2
dev.off()

#WTFD3
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.WTFD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTFD3.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot3 <- FeatureScatter(WTFD3.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(WTFD3.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.WTFD3.tiff", units="in", width=30, height=5, res=300)
plot3 + plot4
dev.off()

#WTFD7
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.WTFD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTFD7.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot5 <- FeatureScatter(WTFD7.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(WTFD7.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.WTFD7.tiff", units="in", width=30, height=5, res=300)
plot5 + plot6
dev.off()

#WTMD0
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.WTMD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTMD0.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(WTMD0.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WTMD0.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.WTM0.tiff", units="in", width=30, height=5, res=300)
plot1 + plot2
dev.off()

#WTMD3
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.WTMD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTMD3.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot3 <- FeatureScatter(WTMD3.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(WTMD3.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.WTMD3.tiff", units="in", width=30, height=5, res=300)
plot3 + plot4
dev.off()

#WTMD7
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot/Vlnplot.WTMD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTMD7.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot5 <- FeatureScatter(WTMD7.SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(WTMD7.SO, feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot/FeatureScatterPlot.WTMD7.tiff", units="in", width=30, height=5, res=300)
plot5 + plot6
dev.off()

###statistical method for determining cutoff
#WT
temp1 <- scater::isOutlier(WTFD0.SO$percent.mt, nmads = 3, type = "higher")
mt.threshold <- min(WTFD0.SO$percent.mt[temp1])
mt.threshold
#  7.321429  - record this value
temp1 <- scater::isOutlier(WTFD3.SO$percent.mt, nmads = 3, type = "higher")
mt.threshold <- min(WTFD3.SO$percent.mt[temp1])
mt.threshold
#  8.979787  - record this value
temp1 <- scater::isOutlier(WTFD7.SO$percent.mt, nmads = 3, type = "higher")
mt.threshold <- min(WTFD7.SO$percent.mt[temp1])
mt.threshold
#  4.125874  - record this value

temp1 <- scater::isOutlier(WTMD0.SO$percent.mt, nmads = 3, type = "higher")
mt.threshold <- min(WTMD0.SO$percent.mt[temp1])
mt.threshold
#  10.12658  - record this value
temp1 <- scater::isOutlier(WTMD3.SO$percent.mt, nmads = 3, type = "higher")
mt.threshold <- min(WTMD3.SO$percent.mt[temp1])
mt.threshold
#  8.830323  - record this value
temp1 <- scater::isOutlier(WTMD7.SO$percent.mt, nmads = 3, type = "higher")
mt.threshold <- min(WTMD7.SO$percent.mt[temp1])
mt.threshold
#  10.79739  - record this value

# Remove these low-quality cells from the dataset using %mt<10
#KO
KOFD0.SO <- subset(KOFD0.SO, subset = nFeature_RNA > 200 & percent.mt < 10)
KOFD3.SO <- subset(KOFD3.SO, subset = nFeature_RNA > 200 & percent.mt < 10)
KOFD7.SO <- subset(KOFD7.SO, subset = nFeature_RNA > 200 & percent.mt < 10)

KOMD0.SO <- subset(KOMD0.SO, subset = nFeature_RNA > 200 & percent.mt < 10)
KOMD3.SO <- subset(KOMD3.SO, subset = nFeature_RNA > 200 & percent.mt < 10)
KOMD7.SO <- subset(KOMD7.SO, subset = nFeature_RNA > 200 & percent.mt < 10)

#WT
WTFD0.SO <- subset(WTFD0.SO, subset = nFeature_RNA > 200 & percent.mt < 10)
WTFD3.SO <- subset(WTFD3.SO, subset = nFeature_RNA > 200 & percent.mt < 10)
WTFD7.SO <- subset(WTFD7.SO, subset = nFeature_RNA > 200 & percent.mt < 10)

WTMD0.SO <- subset(WTMD0.SO, subset = nFeature_RNA > 200 & percent.mt < 10)
WTMD3.SO <- subset(WTMD3.SO, subset = nFeature_RNA > 200 & percent.mt < 10)
WTMD7.SO <- subset(WTMD7.SO, subset = nFeature_RNA > 200 & percent.mt < 10)

# Visualizing the distribution of detected genes, numbers of RNA and the mitochondrial percentage in all cells after removing low-quality cells 
#KOFD0
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.KOFD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOFD0.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(KOFD0.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(KOFD0.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/Vlnplot.KOFD0.tiff", units="in", width=10, height=5, res=300)
plot1 + plot2
dev.off()

#KOFD3
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.KOFD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOFD3.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot3 <- FeatureScatter(KOFD3.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(KOFD3.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/FeatureScatterPlot.KOFD3.tiff", units="in", width=30, height=5, res=300)
plot3 + plot4
dev.off()

#KOFD7
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.KOFD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOFD7.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot5 <- FeatureScatter(KOFD7.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(KOFD7.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/FeatureScatterPlot.KOFD7.tiff", units="in", width=30, height=5, res=300)
plot5 + plot6
dev.off()

#KOMD0
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.KOMD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOMD0.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(KOMD0.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(KOMD0.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/Vlnplot.KOMD0.tiff", units="in", width=10, height=5, res=300)
plot1 + plot2
dev.off()

#KOMD3
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.KOMD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOMD3.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot3 <- FeatureScatter(KOMD3.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(KOMD3.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/FeatureScatterPlot.KOMD3.tiff", units="in", width=30, height=5, res=300)
plot3 + plot4
dev.off()

#KOMD7
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.KOMD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOMD7.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot5 <- FeatureScatter(KOMD7.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(KOMD7.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/FeatureScatterPlot.KOMD7.tiff", units="in", width=30, height=5, res=300)
plot5 + plot6
dev.off()

#WTFD0
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.WTFD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTFD0.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(WTFD0.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WTFD0.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/Vlnplot.WTFD0.tiff", units="in", width=10, height=5, res=300)
plot1 + plot2
dev.off()

#WTFD3
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.WTFD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTFD3.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot3 <- FeatureScatter(WTFD3.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(WTFD3.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/FeatureScatterPlot.WTFD3.tiff", units="in", width=30, height=5, res=300)
plot3 + plot4
dev.off()

#WTFD7
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.WTFD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTFD7.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot5 <- FeatureScatter(WTFD7.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(WTFD7.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/FeatureScatterPlot.WTFD7.tiff", units="in", width=30, height=5, res=300)
plot5 + plot6
dev.off()

#WTMD0
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.WTMD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTMD0.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(WTMD0.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WTMD0.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/Vlnplot.WTMD0.tiff", units="in", width=10, height=5, res=300)
plot1 + plot2
dev.off()

#WTMD3
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.WTMD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTMD3.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot3 <- FeatureScatter(WTMD3.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(WTMD3.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/FeatureScatterPlot.WTMD3.tiff", units="in", width=30, height=5, res=300)
plot3 + plot4
dev.off()

#WTMD7
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Individual Vlnplot_Mt cutoff/Vlnplot.WTMD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTMD7.SO , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot5 <- FeatureScatter(WTMD7.SO , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(WTMD7.SO , feature1 = "nFeature_RNA", feature2 = "percent.mt")
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Ind Feature Scatter Plot_Mt cutoff/FeatureScatterPlot.WTMD7.tiff", units="in", width=30, height=5, res=300)
plot5 + plot6
dev.off()

# Additional important quality control step involves the detection of likely doublets in the dataset. These are cells that were joined during droplet sequencing and will thus result in gene expressions that are not at the single-cell level. 
# For this, we use a tool called scDblFinder. Type in the following commands to run the pipeline:
sce1 <- scDblFinder(GetAssayData(KOFD0.SO , slot="counts"), samples=Idents(KOFD0.SO ))
sce2 <- scDblFinder(GetAssayData(KOFD3.SO , slot="counts"), samples=Idents(KOFD3.SO ))
sce3 <- scDblFinder(GetAssayData(KOFD7.SO , slot="counts"), samples=Idents(KOFD7.SO ))

sce4 <- scDblFinder(GetAssayData(KOMD0.SO , slot="counts"), samples=Idents(KOMD0.SO ))
sce5 <- scDblFinder(GetAssayData(KOMD3.SO , slot="counts"), samples=Idents(KOMD3.SO ))
sce6 <- scDblFinder(GetAssayData(KOMD7.SO , slot="counts"), samples=Idents(KOMD7.SO ))

sce7 <- scDblFinder(GetAssayData(WTFD0.SO , slot="counts"), samples=Idents(WTFD0.SO ))
sce8 <- scDblFinder(GetAssayData(WTFD3.SO , slot="counts"), samples=Idents(WTFD3.SO ))
sce9 <- scDblFinder(GetAssayData(WTFD7.SO , slot="counts"), samples=Idents(WTFD7.SO ))

sce10 <- scDblFinder(GetAssayData(WTMD0.SO , slot="counts"), samples=Idents(WTMD0.SO ))
sce11 <- scDblFinder(GetAssayData(WTMD3.SO , slot="counts"), samples=Idents(WTMD3.SO ))
sce12 <- scDblFinder(GetAssayData(WTMD7.SO , slot="counts"), samples=Idents(WTMD7.SO ))

# Step 2.23: Assign the doublet score to a new metadata variable
KOFD0.SO $scDblFinder.score <- sce1$scDblFinder.score
KOFD3.SO $scDblFinder.score <- sce2$scDblFinder.score
KOFD7.SO $scDblFinder.score <- sce3$scDblFinder.score

KOMD0.SO $scDblFinder.score <- sce4$scDblFinder.score
KOMD3.SO $scDblFinder.score <- sce5$scDblFinder.score
KOMD7.SO $scDblFinder.score <- sce6$scDblFinder.score

WTFD0.SO $scDblFinder.score <- sce7$scDblFinder.score
WTFD3.SO $scDblFinder.score <- sce8$scDblFinder.score
WTFD7.SO $scDblFinder.score <- sce9$scDblFinder.score

WTMD0.SO $scDblFinder.score <- sce10$scDblFinder.score
WTMD3.SO $scDblFinder.score <- sce11$scDblFinder.score
WTMD7.SO $scDblFinder.score <- sce12$scDblFinder.score

# Step 2.24: Visualize the distribution of detected genes, numbers of RNA and the doublet score in all cells
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/KOFD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOFD0.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()

tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/KOFD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOFD3.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()

tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/KOFD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOFD7.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()

tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/KOMD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOMD0.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()

tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/KOMD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOMD3.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()

tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/KOMD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(KOMD7.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()


tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/WTFD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTFD0.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()

tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/WTFD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTFD3.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()

tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/WTFD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTFD7.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()

tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/WTMD0.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTMD0.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()

tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/WTMD3.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTMD3.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()

tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/doublet removal/WTMD7.tiff", units="in", width=10, height=5, res=300)
VlnPlot(WTMD7.SO , features = "scDblFinder.score", raster=FALSE, pt.size=0.5)
dev.off()

# Remove the cells above the 0.25 double score threshold
#KO
WTFD0.SO <- subset(KOFD0.SO , scDblFinder.score < 0.25)
KOFD3.SO <- subset(KOFD3.SO , scDblFinder.score < 0.25)
KOFD7.SO <- subset(KOFD7.SO , scDblFinder.score < 0.25)

KOMD0.SO <- subset(KOMD0.SO , scDblFinder.score < 0.25)
KOMD3.SO <- subset(KOMD3.SO , scDblFinder.score < 0.25)
KOMD7.SO <- subset(KOMD7.SO , scDblFinder.score < 0.25)

#WT
WTFD0.SO <- subset(WTFD0.SO , scDblFinder.score < 0.25)
WTFD3.SO <- subset(WTFD3.SO , scDblFinder.score < 0.25)
WTFD7.SO <- subset(WTFD7.SO , scDblFinder.score < 0.25)

WTMD0.SO <- subset(WTMD0.SO , scDblFinder.score < 0.25)
WTMD3.SO <- subset(WTMD3.SO , scDblFinder.score < 0.25)
WTMD7.SO <- subset(WTMD7.SO , scDblFinder.score < 0.25)
# The dataset is now quality-controlled and it is ready for downstream analysis

# Save the dataset Seurat object as an RDS file into your working directory
saveRDS(KOFD0.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOFD0.SO.post_QC.rds")
saveRDS(KOFD3.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOFD3.SO.post_QC.rds")
saveRDS(KOFD7.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOFD7.SO.post_QC.rds")

saveRDS(KOMD0.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOMD0.SO.post_QC.rds")
saveRDS(KOMD3.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOMD3.SO.post_QC.rds")
saveRDS(KOMD7.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOMD7.SO.post_QC.rds")

saveRDS(WTFD0.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTFD0.SO.post_QC.rds")
saveRDS(WTFD3.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTFD3.SO.post_QC.rds")
saveRDS(WTFD7.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTFD7.SO.post_QC.rds")

saveRDS(WTMD0.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTMD0.SO.post_QC.rds")
saveRDS(WTMD3.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTMD3.SO.post_QC.rds")
saveRDS(WTMD7.SO , file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTMD7.SO.post_QC.rds")

# INTEGRATION: Combining the individual data sets into one post-QC ####
# Single cell data sets are often separated into multiple files because they were sequenced in orig.identityes or groups. In this workflow, we show how to integrate two of the five orig.identityes of the Hu, et al data set. We use the latest methods for dataset integration as described by the following Seurat online vignette (). 
# Make sure all the data sets/samples undergo the same data QC performed in Method 2 (complete - look above)
# Optional Step: If needed, open the  data sets as Seurat objects from their RDS files saved in your working directory
#KOFD0.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOFD0.SO.post_QC.rds")
#KOFD3.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOFD3.SO.post_QC.rds")
#KOFD7.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOFD7.SO.post_QC.rds")

#KOMD0.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOMD0.SO.post_QC.rds")
#KOMD3.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOMD3.SO.post_QC.rds")
#KOMD7.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/KOMD7.SO.post_QC.rds")

#WTFD0.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTFD0.SO.post_QC.rds")
#WTFD3.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTFD3.SO.post_QC.rds")
#WTFD7.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTFD7.SO.post_QC.rds")

#WTMD0.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTMD0.SO.post_QC.rds")
#WTMD3.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTMD3.SO.post_QC.rds")
#WTMD7.SO <-readRDS(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Individual QC/Post-QC Seurat Object/WTMD7.SO.post_QC.rds")

#Adding in additional metadata parameter "bulk identity" and naming 
KOFD0.SO$orig.identity <- KOFD0.SO$bulk.ident
KOFD3.SO$orig.identity <- KOFD3.SO$bulk.ident
KOFD7.SO$orig.identity <- KOFD7.SO$bulk.ident 

KOMD0.SO$orig.identity <- KOMD0.SO$bulk.ident 
KOMD3.SO$orig.identity <- KOMD3.SO$bulk.ident 
KOMD7.SO$orig.identity <- KOMD7.SO$bulk.ident 

WTFD0.SO$orig.identity <- WTFD0.SO$bulk.ident 
WTFD3.SO$orig.identity <- WTFD3.SO$bulk.ident 
WTFD7.SO$orig.identity <- WTFD7.SO$bulk.ident 

WTMD0.SO$orig.identity <- WTMD0.SO$bulk.ident 
WTMD3.SO$orig.identity <- WTMD3.SO$bulk.ident 
WTMD7.SO$orig.identity <- WTMD7.SO$bulk.ident 

KOFD0.SO$bulk.ident <- "KOD0"
KOFD3.SO$bulk.ident <- "KOD3"
KOFD7.SO$bulk.ident <- "KOD7"

KOMD0.SO$bulk.ident <- "KOD0"
KOMD3.SO$bulk.ident <- "KOD3"
KOMD7.SO$bulk.ident <-"KOD7"

WTFD0.SO$bulk.ident <-"WTD0"
WTFD3.SO$bulk.ident <-  "WTD3"
WTFD7.SO$bulk.ident <-  "WTD7"

WTMD0.SO$bulk.ident <- "WTD0"
WTMD3.SO$bulk.ident <-  "WTD3"
WTMD7.SO$bulk.ident <- "WTD7"

#Creating a singular merged file
All.merged<- merge(x = KOFD0.SO , y = c(KOFD3.SO ,KOFD7.SO , KOMD0.SO ,KOMD3.SO , KOMD7.SO ,
                                        WTFD0.SO ,WTFD3.SO ,WTFD7.SO , WTMD0.SO ,WTMD3.SO ,WTMD7.SO ),
                   add.cell.ids = c("KOFD0","KOFD3","KOFD7","KOMD0","KOMD3","KOMD7","WTFD0","WTFD3","WTFD7","WTMD0","WTMD3","WTMD7"), merge.data = TRUE)


# Adding new metadata variable called DPW - D0, D3, D7
All.merged@meta.data$DPW <- "DPW" #make new metadata column

All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "KOFD0")] <- "D0"
All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "KOMD0")] <- "D0"
All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "WTFD0")] <- "D0"
All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "WTMD0")] <- "D0"

All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "KOFD3")] <- "D3"
All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "KOMD3")] <- "D3"
All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "WTFD3")] <- "D3"
All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "WTMD3")] <- "D3"

All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "KOFD7")] <- "D7"
All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "KOMD7")] <- "D7"
All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "WTFD7")] <- "D7"
All.merged@meta.data$DPW[which(All.merged@meta.data$orig.identity == "WTMD7")] <- "D7"

table(All.merged@meta.data$DPW) #checking to see if it was made properly

# Adding new metadata variable "sample.group" - WT vs KO
All.merged@meta.data$sample.group <- "sample.group" #make new metadata column
All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "KOFD0")] <- "KO"
All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "KOMD0")] <- "KO"
All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "WTFD0")] <- "WT"
All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "WTMD0")] <- "WT"

All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "KOFD3")] <- "KO"
All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "KOMD3")] <- "KO"
All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "WTFD3")] <- "WT"
All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "WTMD3")] <- "WT"

All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "KOFD7")] <- "KO"
All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "KOMD7")] <- "KO"
All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "WTFD7")] <- "WT"
All.merged@meta.data$sample.group[which(All.merged@meta.data$orig.identity == "WTMD7")] <- "WT"

table(All.merged@meta.data$sample.group)

# Adding a new metadata variable "Sex": M vs F
All.merged@meta.data$Sex <- "Sex" 

All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "KOFD0")] <- "F"
All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "KOMD0")] <- "M"
All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "WTFD0")] <- "F"
All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "WTMD0")] <- "M"

All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "KOFD3")] <- "F"
All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "KOMD3")] <- "M"
All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "WTFD3")] <- "F"
All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "WTMD3")] <- "M"

All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "KOFD7")] <- "F"
All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "KOMD7")] <- "M"
All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "WTFD7")] <- "F"
All.merged@meta.data$Sex[which(All.merged@meta.data$orig.identity == "WTMD7")] <- "M"

table(All.merged@meta.data$Sex)

# Perform Seurat merging of the  datasets, adding orig.identity-based cell ID annotations 
All.merged<- merge(x = KOFD0.SO , y = c(KOFD3.SO ,KOFD7.SO , KOMD0.SO ,KOMD3.SO , KOMD7.SO ,
                                        WTFD0.SO ,WTFD3.SO ,WTFD7.SO , WTMD0.SO ,WTMD3.SO ,WTMD7.SO ),
                   add.cell.ids = c("KOFD0","KOFD3","KOFD7","KOMD0","KOMD3","KOMD7","WTFD0","WTFD3","WTFD7","WTMD0","WTMD3","WTMD7"), merge.data = TRUE)

#saving the merged file pre-analysis
saveRDS(All.merged, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/All samples merged.rds")

# POST-INTEGRATION: Perform the standard Seurat workflow for the merged datasets: Chose 10 Dims and 0.2 resolution ####
All.merged <- NormalizeData(All.merged)
All.merged <- FindVariableFeatures(All.merged)
All.merged <- ScaleData(All.merged)
All.merged <- RunPCA(All.merged)
ElbowPlot(All.merged, reduction = "pca", ndims = 50) ###will use first 10 PCAs - depicted in dims = 1:10

# Perform UMAP analysis and clustering on the combined dataset prior to data integration
DefaultAssay(All.merged)->"RNA"
All.merged <- RunUMAP(All.merged, dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")
All.merged <- FindNeighbors(All.merged, dims = 1:10, reduction = "pca", graph.name = )
All.merged <- FindClusters(All.merged, resolution = .2, cluster.name = "unintegrated_clusters")

# Visualize the UMAP plot according to cluster and orig.identity numbers 
DimPlot(All.merged, reduction = "umap.unintegrated", group.by = c("seurat_clusters", "orig.identity"))

# Step 7.8: Show distribution of cell numbers in each cluster according to orig.identity number
table(All.merged$orig.identity, All.merged$seurat_clusters)

# Step 7.9: Perform Seurat data integration using the RPCA method --- for more information on this and other data integration methods, please read the following Seurat vignette: 
All.merged <- IntegrateLayers(
  object = All.merged, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE
) ##can do different types integration - just change "method = ..."

# Perform UMAP analysis and clustering on the combined dataset after data integration
All.merged <- RunUMAP(All.merged, reduction = "integrated.rpca", dims = 1:10, reduction.name = "umap.rpca") #choosing 10 dim based on previous
All.merged <- FindNeighbors(All.merged, reduction = "integrated.rpca", dims = 1:10)
All.merged <- FindClusters(All.merged, resolution = .1, cluster.name = "rpca_clusters")

# Visualize the UMAP plot according to cluster and orig.identity numbers after integration: identified 13 clusters - may want to reduce the dimensions to find even less
DimPlot(All.merged, reduction = "umap.rpca", group.by = c("seurat_clusters","orig.identity"))

# Test different resolutions as well. 
#below we are testing different resolutions: 0.3, 0.5, 0.7, and 0.2 ALL WITH SAME NUMBER OF PCAs
All.merged <- FindClusters(All.merged, resolution = .3, cluster.name = "rpca_clusters") ##for iterative studies can do multiple dimensions and resolutions
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/Int.seurat clusters_10 Dim_res 0.3_07122024 PM.tiff", units="in", width=10, height=5, res=300)
DimPlot(All.merged, reduction = "umap.rpca", group.by = c("seurat_clusters","orig.identity"))
dev.off()

All.merged <- FindClusters(All.merged, resolution = .5, cluster.name = "rpca_clusters") ##for iterative studies can do multiple dimensions and resolutions
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/Int.seurat clusters_10 Dim_res 0.5_07122024 PM.tiff", units="in", width=10, height=5, res=300)
DimPlot(All.merged, reduction = "umap.rpca", group.by = c("seurat_clusters","orig.identity"))
dev.off()

All.merged <- FindClusters(All.merged, resolution = .7, cluster.name = "rpca_clusters") ##for iterative studies can do multiple dimensions and resolutions
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/Int.seurat clusters_10 Dim_res 0.7_07162024 AM.tiff", units="in", width=10, height=5, res=300)
DimPlot(All.merged, reduction = "umap.rpca", group.by = c("seurat_clusters","orig.identity"))
dev.off()

All.merged <- FindClusters(All.merged, resolution = .2, cluster.name = "rpca_clusters") ##for iterative studies can do multiple dimensions and resolutions
tiff("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/Int.seurat clusters_10 Dim_res 0.2_07122024 PM.tiff", units="in", width=10, height=5, res=300)
DimPlot(All.merged, reduction = "umap.rpca", group.by = c("seurat_clusters","orig.identity"))
dev.off()

#will choose 10 dimensions with resolution =0.2
table(All.merged$orig.identity, All.merged$seurat_clusters)
All.merged <- JoinLayers(All.merged)

# From the UMAP plot and the table of the integrated data, there is now excellent overlap between the two orig.identityes across different clusters. Interestingly, after integration the numbers of clusters at this resolution is the same as what was determined for orig.identity #1 alone.
# After dataset integration and prior to downstream analyses, the layers of the merged dataset must be joined
All.merged <- JoinLayers(All.merged)

# Save the dataset Seurat object as an RDS file into your working directory
saveRDS(All.merged, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/All.merged.postM3-81024.rds")

# POST-INTEGRATION: Trying other PCAs and resolutions to find best cluster identification (post-integration) + Marker GEne Identification ####
# Note: if rerunning from an already integrated dataset (like I am...use this code below)
DefaultAssay(All.merged)<-"SCT"  #let's R know we are using integrated data set
All.merged <- NormalizeData(All.merged)
All.merged <- FindVariableFeatures(All.merged)
All.merged <- ScaleData(All.merged)
All.merged <- RunPCA(All.merged)
All.merged <- RunUMAP(All.merged, reduction = "integrated.rpca", dims = 1:10, reduction.name = "umap.rpca") #choosing 10 dim based on previous
All.merged <- FindNeighbors(All.merged, reduction = "integrated.rpca", dims = 1:10)
All.merged <- FindClusters(All.merged, resolution = .2, cluster.name = "rpca_clusters", graph.name = "RNA_snn") #need to provide graphname

All.merged <- RunUMAP(All.merged, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca") #choosing 10 dim based on previous
All.merged <- FindNeighbors(All.merged, reduction = "integrated.rpca", dims = 1:15)
All.merged <- FindClusters(All.merged, resolution = .4, cluster.name = "rpca_clusters", graph.name = "RNA_snn")

# Visualize the UMAP plot according to cluster and orig.identity numbers after integration: identified 13 clusters - may want to reduce the dimensions to find even less
DimPlot(All.merged, reduction = "umap.rpca", group.by = c("seurat_clusters"), label = TRUE)
#DimPlot(All.merged, reduction = "umap.rpca", group.by = c("seurat_clusters","orig.identity"), label = TRUE)

# Show distribution of cell numbers in each cluster according to orig.identity number after integration
table(All.merged$orig.identity, All.merged$seurat_clusters)

# From the UMAP plot and the table of the integrated data, there is now excellent overlap between the two orig.identityes across different clusters. Interestingly, after integration the numbers of clusters at this resolution is the same as what was determined for orig.identity #1 alone.
# After dataset integration and prior to downstream analyses, the layers of the merged dataset must be joined
All.merged <- JoinLayers(All.merged, assay = "RNA")

# Save the dataset Seurat object as an RDS file into your working directory
saveRDS(All.merged, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/All.merged.postM3-81024.rds") #15 PCA Dimensions 0.4 Resolutions
saveRDS(All.merged, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/All.Merged.post method 3_08272024.rds") #10 PCA Dimensions 0.2 Resolutions

#The RDS file I was using for earlier analysis (the one shown in most recent thesis update is "All.Merged.post method 3_08062024.rds") 
#This one was with 10dim0.2res 
# Optional: Reloading saved RDS
All.merged <- readRDS("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/All.merged.postM3-81024.rds")
All.merged <- readRDS("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/All.Merged.post method 3_08062024.rds")

# Determine identities of the cell clusters
Idents(All.merged) <- "seurat_clusters"
Cell_markers.All <- FindAllMarkers(All.merged, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
write.csv(Cell_markers.All, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/15Dim-0.4Res.MarkerGenes.81024.txt") 
#here we are utilizing 10dim_0.2res for seurat clusters (so 17 clusters)

# Assigning cell identities for 10 dim and 0.2 res: based on EnrichR of top genes ####
Idents(All.merged) <- "seurat_clusters"
All.merged[["cell_types"]] <- Idents(All.merged)
Idents(All.merged) <- "cell_types"
All.merged <- RenameIdents(All.merged, 
                           "0" = "Macrophages",	
                           "1" = "Fibroblasts",	
                           "2" = "Neutrophils",	
                           "3" = "T Cells",	
                           "4" = "Keratinocytes",	
                           "5" = "Keratinocytes",	
                           "6" = "Dendritic Cells",	
                           "7" = "T Cells",	
                           "8" = "Keratinocytes",	
                           "9" = "Keratinocytes",	
                           "10" = "Fibroblasts",	
                           "11" = "Endothelial Cells",	
                           "12" = "T Cells",	
                           "13" = "Smooth Muscle Cells",	
                           "14" = "Fibroblasts",	
                           "15" = "Pericytes",	
                           "16" = "Endothelial Cells",	
                           "17" = "Keratinocytes"	
)
levels(All.merged)
All.merged[["cell_types"]] <- Idents(All.merged)

# Visualize the clusters with cell types
DimPlot(All.merged, reduction = "umap.rpca", group.by = c("cell_types"), label = TRUE, repel = TRUE)

#Potential Markers for below
#C3ar1 Macrophage
#Fstl1, Ncam1, Perp, Rarres2, Col6a1, Prkg1, FB
#S100a9 Neutrophil
#Cd3e, Skap1, Stmn1, T-cells
#Dsp, Krtdap, Krt14, Epithelial Cells
#Cd83 Dendritic Cells 
#Pecam1, Gng11 Endothelial Cells
#Apoc1, ldhb Adipocytes


# Visualize the localization of the top (bolded) cluster marker genes used to characterize each major cell type 
# this will have multiple feature plots - one for each cell group with a cluster of markers 
# Fibroblast
FB_markers_Mfap5 <- list(c("Col5a1","Fstl1","Lum","Dcn","Mmp2"))
All.merged <- AddModuleScore(object = All.merged, features = FB_markers_Mfap5, name = "FB Markers")
FeaturePlot(All.merged, features = "FB Markers1",raster=FALSE, reduction = "umap.rpca")
# FeaturePlot(All.merged, features = c("Col5a1","Fstl1","Lum","Dcn","Mmp2"),raster=FALSE, reduction = "umap.rpca")

# Macrophage
Macrophage_markers_Mfap5 <- list(c("Ctss","Pid1","Mpeg1","Ccl9","Cd68"))
All.merged <- AddModuleScore(object = All.merged, features = Macrophage_markers_Mfap5, name = "Macrophage Markers")
FeaturePlot(All.merged, features = "Macrophage Markers1",raster=FALSE, reduction = "umap.rpca")
# FeaturePlot(All.merged, features = c("Ctss","Pid1","Mpeg1","CCl9","Cd68"),raster=FALSE, reduction = "umap.rpca")

# Neutrophil
Neutrophil_markers_Mfap5 <- list(c("S100a9","S100a8","Slc7a11","Il1b","Hdc"))
All.merged <- AddModuleScore(object = All.merged, features = Neutrophil_markers_Mfap5, name = "Neutrophil Markers")
FeaturePlot(All.merged, features = "Neutrophil Markers1",raster=FALSE, reduction = "umap.rpca")
# FeaturePlot(All.merged, features = ("S100a9","S100a8","Slc7a11","Il1b","Hdc"),raster=FALSE, reduction = "umap.rpca")

# T cell
Tcell_markers_Mfap5 <- list(c("Cd3e","Cd247","Tox","Skap1","Ikzf3"))
All.merged <- AddModuleScore(object = All.merged, features = Tcell_markers_Mfap5, name = "T Cell Markers")
FeaturePlot(All.merged, features = "T Cell Markers1",raster=FALSE, reduction = "umap.rpca")
# FeaturePlot(All.merged, features = "c("Cd3e","Cd247","Tox","Skap1","Ikzf3")",raster=FALSE, reduction = "umap.rpca")

# Endothelial cell
EndothelialCell_markers_Mfap5 <- list(c("Pecam1","Cdh5","Col4a1","Gng11","Cdh5"))
All.merged <- AddModuleScore(object = All.merged, features = EndothelialCell_markers_Mfap5, name = "Endothelial Cell Markers")
FeaturePlot(All.merged, features = "Endothelial Cell Markers1",raster=FALSE, reduction = "umap.rpca")
# FeaturePlot(All.merged, features = "c("Pecam1","Cdh5","Col4a1","Gng11","Cdh5")",raster=FALSE, reduction = "umap.rpca")

# Smooth Muscle cell
SmoothMuscle_marker_Mfap5 <- list(c("Ncam1","Pax7","Drp2","Megf10","Cdh15"))
All.merged <- AddModuleScore(object = All.merged, features = SmoothMuscle_marker_Mfap5, name = "SmoothMuscle Markers")
FeaturePlot(All.merged, features = "SmoothMuscle Markers1",raster=FALSE, reduction = "umap.rpca")
# FeaturePlot(All.merged, features = "c("Ncam1","Pax7","Drp2","Megf10","Cdh15")",raster=FALSE, reduction = "umap.rpca")

# Pericyte
Pericyte_marker_Mfap5 <- list(c("Prkg1","Mylk","Acta2","Sparcl1","Notch3"))
All.merged <- AddModuleScore(object = All.merged, features = Pericyte_marker_Mfap5, name = "Pericyte Markers")
FeaturePlot(All.merged, features = "Pericyte Markers1",raster=FALSE, reduction = "umap.rpca")
# FeaturePlot(All.merged, features = "c("Prkg1","Mylk","Acta2","Sparcl1","Notch3")",raster=FALSE, reduction = "umap.rpca")

# Keratinocyte
Keratinocyte_marker_Mfap5 <- list(c("Krtdap","Dsp","Apoc1","Krt1","Krt5"))
All.merged <- AddModuleScore(object = All.merged, features = Keratinocyte_marker_Mfap5, name = "Keratinocyte Markers")
FeaturePlot(All.merged, features = "Keratinocyte Markers1",raster=FALSE, reduction = "umap.rpca")
# FeaturePlot(All.merged, features = "c("Krtdap","Dsp","Apoc1","Krt1","Krt5")",raster=FALSE, reduction = "umap.rpca")

# Dendritic Cell
Dendritic_marker_Mfap5 <- list(c("Cd83","Flt3","P2ry10","Tbc1d4","Rtn1"))
All.merged <- AddModuleScore(object = All.merged, features = Dendritic_marker_Mfap5, name = "Dendritic Markers")
FeaturePlot(All.merged, features = "Dendritic Markers1",raster=FALSE, reduction = "umap.rpca")
# FeaturePlot(All.merged, features = "c("Cd83","Flt3","P2ry10","Tbc1d4","Rtn1")",raster=FALSE, reduction = "umap.rpca")

# Perform SCTransform which transforms the gene expression values to better visualize DEGs on a dot plot
future::plan("multisession", workers = 5) # do parallel
options(future.globals.maxSize = 48000 * 1024^2) #increases the max global export for future expression 

DefaultAssay(All.merged) <- "RNA"
All.merged <- SCTransform(object = All.merged, assay = "RNA")

#10 PCA Dimensions 0.2 Resolutions : 18 clusters
#labeled the cell types
#LOAD THIS IF YOU NEED TO RELOAD SO YOU CAN DO ANALYSES ON THIS SEURAT OBJECT...I AM USING THIS
#All.merged <- readRDS("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/All.merged.postM3-82724.rds")

cell.proportion<-table(All.merged$bulk.ident, All.merged$cell_types)
write.csv(cell.proportion, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/CellProportions.txt")

#saveRDS(All.merged, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/All.Merged.post method 3_08272024") 
#15 PCA Dimensions 0.4 Resolutions : 28 clusters

# Visualize the localization of the top (bolded) cluster marker genes on a series of UMAP plots - all seurat clusters (18 UMAPS)
DefaultAssay(All.merged) <- "RNA"
FeaturePlot(All.merged, features = c("Fstl1","Ctss","S100a9","Cd3e","Pecam1","Ncam1","Acta2",
                                     "Krtdap","Cd83"), raster=FALSE, ncol = 5, reduction = "umap.rpca")


# Visualize the top cluster marker DEGs on a dot plot, grouped by the original cluster numbers and the annotated cell types
DotPlot(All.merged, group.by = "seurat_clusters", features = c("Fstl1","Ctss","S100a9","Cd3e","Pecam1","Ncam1","Acta2",
                                                               "Krtdap","Cd83")) + RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

DotPlot(All.merged, group.by = "cell_types", features = c("Fstl1","Ctss","S100a9","Cd3e","Pecam1","Ncam1","Acta2",
                                                          "Krtdap","Cd83")) + RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

# Two dot plots confirming the high level of expression of the top cell marker genes only in their respective seurat clusters and major cell types

#Looking at the change in cell type proportion over time 
pt1 <- table(All.merged$DPW, All.merged$cell_types)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)
ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") + RotatedAxis() +
  ylab("Proportion") + 
  scale_fill_manual(values = brewer.pal(8, "Paired")) +
  theme(legend.title = element_blank())

pt2 <- table(All.merged$cell_types, All.merged$DPW)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") + 
  ylab("Proportion") + 
  scale_fill_manual(values = brewer.pal(9, "Paired")) +
  theme(legend.title = element_blank())

pt3 <- table(All.merged$cell_types, All.merged$bulk.ident)
pt3<-pt3[,c(1,4,2,5,3,6)]
write.xlsx(pt3, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/Table of Cell Type and Frequency.xlsx")
pt3 <- as.data.frame(pt3)
pt3$Var1 <- as.character(pt3$Var1)
ggplot(pt3, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") + 
  ylab("Proportion") + 
  scale_fill_manual(values = brewer.pal(9, "Paired")) +
  theme(legend.title = element_blank())

pt4 <- table(All.merged$cell_types, All.merged$orig.identity)
pt4 <- as.data.frame(pt4)
pt4$Var1 <- as.character(pt4$Var1)
ggplot(pt4, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") + 
  ylab("Proportion") + 
  scale_fill_manual(values = brewer.pal(9, "Paired")) +
  theme(legend.title = element_blank())

table(All.merged$cell_types, All.merged$orig.identity)

saveRDS(All.merged, "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/MFAP5 Paper 2/GEO submission_scRNA seq files/Processed Data Files/Merged Dataset_05122025.rds")

# FB subpopulation identification and analyses  ####
# Optional step: If needed, load in the saved RDS file as a Seurat object
# All.merged <- readRDS("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/All.merged.postM3-82724")
# Subset the original dataset according to the fibroblast cell identity
Idents(All.merged) <- "cell_types"
dataset_fibroblast <- subset(All.merged, idents = "Fibroblasts")

# Perform PCA on this smaller dataset and visualize the amount of dataset variation with respect to PCA dimensions
dataset_fibroblast <- RunPCA(dataset_fibroblast, verbose = TRUE)
ElbowPlot(dataset_fibroblast, reduction = "pca", ndims = 50)

# Perform UMAP dimensional reduction and finding neighbors analysis using the first 5 PCA dimensions
dataset_fibroblast <- RunUMAP(dataset_fibroblast, verbose = TRUE, dims = 1:15)
dataset_fibroblast <- FindNeighbors(dataset_fibroblast, verbose = TRUE, dims = 1:15)

# Perform cell clustering of the dataset using a set conservative resolution of 0.1, but check others
dataset_fibroblast <- FindClusters(dataset_fibroblast, verbose = TRUE, resolution = 0.1)

# Visualize the clustering of the cells on a UMAP plot
DimPlot(dataset_fibroblast, reduction = "umap.rpca", group.by = c("seurat_clusters"), raster = FALSE, label = TRUE, pt.size = 0.5, label.size =5)

# using resolution of 0.1 we get 6 clusters (0-5)
# For sake of completeness, will check other resolutions
# dataset_fibroblast <- FindClusters(dataset_fibroblast, verbose = TRUE, resolution = 0.2)
# DimPlot(dataset_fibroblast, reduction = "umap.rpca", group.by = c("seurat_clusters"), raster = FALSE, label = TRUE, pt.size = 0.5, label.size =5) 
    #9 clusters...with 0 and 1 being very similar...so perhaps over clustering?

# dataset_fibroblast <- FindClusters(dataset_fibroblast, verbose = TRUE, resolution = 0.4)
# DimPlot(dataset_fibroblast, reduction = "umap.rpca", group.by = c("seurat_clusters"), raster = FALSE, label = TRUE, pt.size = 0.5, label.size =5)
    #12 clusters...definitely appears to have some over clustering 

# Save the dataset Seurat object as an RDS file into your working directory
saveRDS(dataset_fibroblast, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/FBdataset_0-5 clusters_090924.rds")

# Visualize the cells on a UMAP plot
DimPlot(dataset_fibroblast, reduction = "umap.rpca", group.by = "bulk.ident", raster = FALSE, label = TRUE, pt.size = 0.5, label.size = 3, repel = TRUE)
DimPlot(dataset_fibroblast, reduction = "umap.rpca", split.by = "bulk.ident", raster = FALSE, label = TRUE, pt.size = 0.25, label.size = 3, repel = TRUE)
DimPlot(dataset_fibroblast, reduction = "umap.rpca", split.by = "DPW", raster = FALSE, label = TRUE, pt.size = 0.25, label.size = 3, repel = TRUE)
DimPlot(dataset_fibroblast, reduction = "umap.rpca", group.by = "bulk.ident",split.by = "RNA_snn_res.0.1" , raster = FALSE, label = TRUE, pt.size = 0.25, label.size = 3, repel = TRUE)

# Visualize the wound timecourse annotation of the cells on a UMAP plot based on sex
DimPlot(dataset_fibroblast, reduction = "umap.rpca", group.by = "Sex", split.by = "bulk.ident", raster = FALSE, label = TRUE, pt.size = 0.25, label.size = 3, repel = TRUE)
DimPlot(dataset_fibroblast, reduction = "umap.rpca", group.by = "Sex", split.by = "RNA_snn_res.0.1", raster = FALSE, label = TRUE, pt.size = 0.25, label.size = 3, repel = TRUE)

# Generate tables of how many cells of each type occur in each DPW
table(dataset_fibroblast$seurat_clusters, dataset_fibroblast$DPW)

# Visualize the proportion of DPW in each cell type
pt1 <- table(dataset_fibroblast$seurat_clusters, dataset_fibroblast$DPW)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)
ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") + RotatedAxis() +
  ylab("Proportion") + 
  scale_fill_manual(values = brewer.pal(8, "Paired")) +
  theme(legend.title = element_blank())

pt1 <- table(dataset_fibroblast$DPW, dataset_fibroblast$ExpGroup)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)
ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") + RotatedAxis() +
  ylab("Proportion") + 
  scale_fill_manual(values = brewer.pal(8, "Paired")) +
  theme(legend.title = element_blank())

pt2 <-  (table(dataset_fibroblast$seurat_clusters, dataset_fibroblast$bulk.ident))
pt2<-pt2[,c(1,4,2,5,3,6)]
write.xlsx(pt2, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Table of FB Frequency.xlsx")
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") + 
  ylab("Proportion") + 
  scale_fill_manual(values = brewer.pal(9, "Paired")) +
  theme(legend.title = element_blank())

pt2 <-  table(dataset_fibroblast$seurat_clusters, dataset_fibroblast$orig.identity)
pt2<-pt2[,c(1,4,7,10,2,5,8,11,3,6,9,12)]
write.xlsx(pt2, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Table of FB Frequency by ExpGroup.xlsx")
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") + 
  ylab("Proportion") + 
  scale_fill_manual(values = brewer.pal(9, "Paired")) +
  theme(legend.title = element_blank())

saveRDS(dataset_fibroblast, "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/MFAP5 Paper 2/GEO submission_scRNA seq files/Processed Data Files/FB subclusters_0-6_05122025.rds")

# Get DEG lists for the 6 fibroblast subtypes and save them into a text file in the working directory
Idents(dataset_fibroblast) <- "seurat_clusters"
fibroblast_markers <- FindAllMarkers(dataset_fibroblast, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(fibroblast_markers, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/090624_dim15res0.1_FB_markers.txt")
#this will get you marker genes for each FB group so you can make feature plots with module scoring of specific FB marker genes
#this is helpful for identifying a FB subpopulation that might have MFAP5 has a marker gene

# Creating module scoring feature plots for the 6 FB subclusters
# FB Cluster 0
FB_markers_C0<- list(c("Cthrc1","Tnc","Lrrc15","Col12a1","Mfap5"))
dataset_fibroblast <- AddModuleScore(object = dataset_fibroblast, features = FB_markers_C0, name = "FB C0")
FeaturePlot(dataset_fibroblast, features = "FB C01",raster=FALSE, reduction = "umap.rpca")

# FB Cluster 1
FB_markers_C1<- list(c("Scara5","Entpd2","Ltbp4","Cadm3","Deptor"))
dataset_fibroblast <- AddModuleScore(object = dataset_fibroblast, features = FB_markers_C1, name = "FB C1")
FeaturePlot(dataset_fibroblast, features = "FB C11",raster=FALSE, reduction = "umap.rpca")

# FB Cluster 2
FB_markers_C2<- list(c("Smpd3","Scara5","Akr1c18","Dpp4","Mfap5"))
dataset_fibroblast <- AddModuleScore(object = dataset_fibroblast, features = FB_markers_C2, name = "FB C2")
FeaturePlot(dataset_fibroblast, features = "FB C21",raster=FALSE, reduction = "umap.rpca")

# FB Cluster 3
FB_markers_C3<- list(c("Timp1","Hif1a","Tagln2"))
dataset_fibroblast <- AddModuleScore(object = dataset_fibroblast, features = FB_markers_C3, name = "FB C3")
FeaturePlot(dataset_fibroblast, features = "FB C31",raster=FALSE, reduction = "umap.rpca")

# FB Cluster 4
FB_markers_C4<- list(c("Vit","Tenm2","Enpp2","Vwa1"))
dataset_fibroblast <- AddModuleScore(object = dataset_fibroblast, features = FB_markers_C4, name = "FB C4")
FeaturePlot(dataset_fibroblast, features = "FB C41",raster=FALSE, reduction = "umap.rpca")

# FB Cluster 5
FB_markers_C5<- list(c("Acta1","Myl1","Tnnc2","Tnni2","Ckm"))
dataset_fibroblast <- AddModuleScore(object = dataset_fibroblast, features = FB_markers_C5, name = "FB C5")
FeaturePlot(dataset_fibroblast, features = "FB C51",raster=FALSE, reduction = "umap.rpca")

# Visualize the genes in the list in the fibroblast-only dataset by calling the variable in the features parameter of the dotplot
DefaultAssay(dataset_fibroblast) <- "SCT"
FB_type_marker = c("Cthrc1","Tnc","Lrrc15","Col12a1","Mfap5","Scara5","Entpd2","Ltbp4","Cadm3","Deptor",
                   "Smpd3","Akr1c18","Dpp4","Mylk","Timp1","Hif1a","Tagln2","Vit","Tenm2","Enpp2","Vwa1",
                   "Acta1","Myl1","Tnnc2","Tnni2","Ckm" )
DotPlot(dataset_fibroblast, group.by="seurat_clusters", features = FB_type_marker) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "grey", high = "red")

# Visualize the genes in the list in the original single-cell dataset by calling the variable in the features parameter of the dotplot
DefaultAssay(All.merged) <- "SCT"
FB.only.markers = c("Cthrc1","Tnc","Lrrc15","Col12a1","Mfap5","Sfrp2","Col1a1","Col3a1")
DotPlot(All.merged, group.by="cell_types", features = FB.only.markers) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "grey", high = "red")
DotPlot(All.merged, group.by="cell_types", features = FB_type_marker) + RotatedAxis() + scale_colour_gradient2(midpoint = 0, low = "blue", mid = "grey", high = "red")

#Additional DimPlots
DimPlot(dataset_fibroblast, reduction = "umap.rpca", group.by = "Group", split.by = "DPW", label = TRUE, label.size = 2.5, pt.size = 0.25, repel = TRUE)
DimPlot(dataset_fibroblast, reduction = "umap.rpca", group.by = "DPW", split.by = "Group", label = TRUE, label.size = 2.5, pt.size = 0.25, repel = TRUE)
DimPlot(dataset_fibroblast, reduction = "umap.rpca", group.by = c("DPW", "Group"), label = TRUE, label.size = 2.5, pt.size = 0.25, repel = TRUE)

# Differential Gene Expression comparing WT and KO (non-pseudobulk)
table(dataset_fibroblast$ExpGroup)
Idents(dataset_fibroblast) <- "ExpGroup"
FBs_markers.WTvsKO <- FindAllMarkers(dataset_fibroblast, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.DEGs.txt")

# Differential Gene Expression comparing WT and KO (non-pseudobulk): Stratified across time points
Idents(dataset_fibroblast) <- "DPW"
table(dataset_fibroblast.D0$ExpGroup)
dataset_fibroblast.D0<-subset(dataset_fibroblast, ident = "D0")
dataset_fibroblast.D3<-subset(dataset_fibroblast, ident = "D3")
dataset_fibroblast.D7<-subset(dataset_fibroblast, ident = "D7")

Idents(dataset_fibroblast.D0) <- "ExpGroup"
FBs_markers.WTvsKO.D0 <- FindAllMarkers(dataset_fibroblast.D0, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.D0, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.D0.txt")

Idents(dataset_fibroblast.D3) <- "ExpGroup"
FBs_markers.WTvsKO.D3 <- FindAllMarkers(dataset_fibroblast.D3, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.D3, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.D3.txt")

Idents(dataset_fibroblast.D7) <- "ExpGroup"
FBs_markers.WTvsKO.D7 <- FindAllMarkers(dataset_fibroblast.D7, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.D7, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.D7.txt")

# Labeling the FBs into subgroups : Creating a metadata variable could help with future pseudobulk analysis ####
#If we want to start to compare WT and KO FBs within the same group, we will have to split it
dataset_fibroblast@meta.data$FBgroup <-"FBgroup" #make new metadata column - groups WT and KO together based on DPW
#Labeling KOD0 FBs into subgroups 0-5
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD0" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "KOD0-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD0" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "KOD0-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD0" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "KOD0-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD0" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "KOD0-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD0" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "KOD0-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD0" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "KOD0-5"

dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD0" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "KOD0-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD0" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "KOD0-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD0" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "KOD0-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD0" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "KOD0-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD0" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "KOD0-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD0" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "KOD0-5"

#Labeling KOD3 FBs into subgroups 0-5
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD3" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "KOD3-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD3" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "KOD3-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD3" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "KOD3-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD3" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "KOD3-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD3" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "KOD3-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD3" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "KOD3-5"

dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD3" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "KOD3-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD3" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "KOD3-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD3" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "KOD3-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD3" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "KOD3-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD3" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "KOD3-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD3" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "KOD3-5"

#Labeling KOD7 FBs into subgroups 0-5
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD7" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "KOD7-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD7" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "KOD7-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD7" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "KOD7-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD7" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "KOD7-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD7" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "KOD7-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOFD7" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "KOD7-5"

dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD7" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "KOD7-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD7" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "KOD7-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD7" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "KOD7-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD7" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "KOD7-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD7" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "KOD7-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "KOMD7" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "KOD7-5"

#Labeling WTD0 FBs into subgroups 0-5
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD0" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "WTD0-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD0" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "WTD0-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD0" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "WTD0-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD0" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "WTD0-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD0" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "WTD0-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD0" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "WTD0-5"

dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD0" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "WTD0-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD0" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "WTD0-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD0" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "WTD0-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD0" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "WTD0-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD0" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "WTD0-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD0" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "WTD0-5"

#Labeling WTD3 FBs into subgroups 0-5
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD3" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "WTD3-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD3" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "WTD3-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD3" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "WTD3-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD3" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "WTD3-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD3" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "WTD3-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD3" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "WTD3-5"

dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD3" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "WTD3-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD3" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "WTD3-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD3" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "WTD3-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD3" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "WTD3-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD3" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "WTD3-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD3" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "WTD3-5"

#Labeling WTD7 FBs into subgroups 0-5
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD7" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "WTD7-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD7" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "WTD7-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD7" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "WTD7-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD7" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "WTD7-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD7" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "WTD7-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTFD7" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "WTD7-5"

dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD7" & dataset_fibroblast@meta.data$seurat_clusters == "0")] <- "WTD7-0"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD7" & dataset_fibroblast@meta.data$seurat_clusters == "1")] <- "WTD7-1"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD7" & dataset_fibroblast@meta.data$seurat_clusters == "2")] <- "WTD7-2"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD7" & dataset_fibroblast@meta.data$seurat_clusters == "3")] <- "WTD7-3"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD7" & dataset_fibroblast@meta.data$seurat_clusters == "4")] <- "WTD7-4"
dataset_fibroblast@meta.data$FBgroup[which(dataset_fibroblast@meta.data$orig.identity == "WTMD7" & dataset_fibroblast@meta.data$seurat_clusters == "5")] <- "WTD7-5"

table(dataset_fibroblast$FBgroup)
# Fibroblasts Subpopulation identification and Analysis for MFAP5+ expressing FBs  ####
# If needed, load in the saved RDS file as a Seurat object
# All.merged<- readRDS("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/All.merged.postM3-82724.rds")

# Subset the original dataset according to the fibroblast cell identity
Idents(All.merged) <- "cell_types"
table(All.merged$cell_types, All.merged$bulk.ident) #less smooth muscles in KO compared to WT at D0 (about equal throughout healing)

ds_FBs <- subset(All.merged, idents = "Fibroblasts")

ds_FBs@meta.data$ExpGroup<-"ExpGroup"
ds_FBs@meta.data$ExpGroup[which(ds_FBs@meta.data$bulk.ident == "KOD0")] <- "KO"
ds_FBs@meta.data$ExpGroup[which(ds_FBs@meta.data$bulk.ident == "KOD3")] <- "KO"
ds_FBs@meta.data$ExpGroup[which(ds_FBs@meta.data$bulk.ident == "KOD7")] <- "KO"
ds_FBs@meta.data$ExpGroup[which(ds_FBs@meta.data$bulk.ident == "WTD0")] <- "WT"
ds_FBs@meta.data$ExpGroup[which(ds_FBs@meta.data$bulk.ident == "WTD3")] <- "WT"
ds_FBs@meta.data$ExpGroup[which(ds_FBs@meta.data$bulk.ident == "WTD7")] <- "WT"
table(ds_FBs$ExpGroup)

saveRDS(
  object = ds_FBs,
  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/RDS Files/Fibroblasts Dataset.rds")

# Perform PCA on this smaller dataset and visualize the amount of dataset variation with respect to PCA dimensions
ds_FBs <- RunPCA(ds_FBs, verbose = TRUE)
ElbowPlot(ds_FBs, reduction = "pca", ndims = 50)
# Much of the major variation occurs within the first 10 dimensions
# Perform UMAP dimensional reduction and finding neighbors analysis using the first 7 PCA dimensions
ds_FBs <- RunUMAP(ds_FBs, verbose = TRUE, dims = 1:15)
ds_FBs <- FindNeighbors(ds_FBs, verbose = TRUE, dims = 1:15)
#When I performed using Dim = 10, there were 5 seurat clusters, but cluster 5 had so few cells so changed to 7

# Perform cell clustering of the dataset using a set conservative resolution of 0.1, but check others
ds_FBs <- FindClusters(ds_FBs, verbose = TRUE, resolution = 0.1) # using resolution of 0.1 we get 6 clusters (0-6)

# Visualize the clustering of the cells on a UMAP plot
DimPlot(ds_FBs, reduction = "umap.rpca", group.by = c("seurat_clusters"), raster = FALSE, label = TRUE, pt.size = 1, label.size =5)
DimPlot(ds_FBs, reduction = "umap.rpca", group.by = c("seurat_clusters"), 
        split.by = "bulk.ident", raster = FALSE, label = TRUE, pt.size = 1.5, label.size =5)
DimPlot(ds_FBs, reduction = "umap.rpca", group.by = "ExpGroup", raster = FALSE, label = TRUE, pt.size = 1, label.size =5, repel=TRUE)
DimPlot(ds_FBs, reduction = "umap.rpca", group.by = "bulk.ident", raster = FALSE, label = TRUE, pt.size = 1, label.size = 3, repel = TRUE)
DimPlot(ds_FBs, reduction = "umap.rpca", group.by = "ExpGroup", split.by = "DPW", raster = FALSE, label = TRUE, pt.size = 1, label.size = 3, repel = TRUE)

# Visualize the wound timecourse annotation of the cells on a UMAP plot based on sex
DimPlot(ds_FBs, reduction = "umap.rpca", group.by = "Sex", split.by = "bulk.ident", raster = FALSE, label = TRUE, pt.size = 1, label.size = 3, repel = TRUE)

# Generate tables of how many cells of each type occur in each DPW
table(ds_FBs$DPW, ds_FBs$ExpGroup)
table(ds_FBs$seurat_clusters, ds_FBs$bulk.ident)

# Visualize the proportion of DPW in each cell type
pt1 <- table(ds_FBs$DPW, ds_FBs$ExpGroup)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt1$Var1)
ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") + RotatedAxis() +
  ylab("Proportion") + 
  scale_fill_manual(values = brewer.pal(8, "Paired")) +
  theme(legend.title = element_blank())

pt2 <-  (table(ds_FBs$seurat_clusters, ds_FBs$bulk.ident))
pt2<-pt2[,c(1,4,2,5,3,6)]
write.xlsx(pt2, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/SMC Analysis/Table of SMC Frequency.xlsx")
pt2 <- as.data.frame(pt2)

pt2$Var1 <- as.character(pt2$Var1)
ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") + 
  ylab("Proportion") + 
  scale_fill_manual(values = brewer.pal(9, "Paired")) +
  theme(legend.title = element_blank())

pt2 <-  table(ds_FBs$seurat_clusters, ds_FBs$orig.identity)
pt2<-pt2[,c(1,4,7,10,2,5,8,11,3,6,9,12)]
write.xlsx(pt2, "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/SMC Analysis/Table of SMC Frequency by ExpGroup.xlsx")
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") + 
  ylab("Proportion") + 
  scale_fill_manual(values = brewer.pal(9, "Paired")) +
  theme(legend.title = element_blank())


# Get DEG lists for the 6 FB subtypes and save them into a text file in the working directory
Idents(ds_FBs) <- "seurat_clusters"
FBs_markers <- FindAllMarkers(ds_FBs, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_FB_markers.txt")
#fibroblast subclusters 0 and 2 have MFAP5 as a marker gene

# EnrichR analysis + Generation of bubble plots for clusters 0 and 2
#EnrichR Analysis
FBs_markers.C0 <- FBs_markers %>%filter(p_val_adj <  0.01 & avg_log2FC > 1 & pct.1 > 0.4 & cluster == "0")
FBs_markers.C2 <- FBs_markers %>%filter(p_val_adj <  0.01 & avg_log2FC > 1 & pct.1 > 0.4 & cluster == "2")

EnrichR.FBs.C0<- enrichr(FBs_markers.C0$gene, dbs)
EnrichR.FBs.C2<- enrichr(FBs_markers.C2$gene, dbs)

#saving enrichment analysis 
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP C0")
addWorksheet(Enrichments, "MF C0")
addWorksheet(Enrichments, "CC C0")
addWorksheet(Enrichments, "Reactome C0")
addWorksheet(Enrichments, "BP C2")
addWorksheet(Enrichments, "MF C2")
addWorksheet(Enrichments, "CC C2")
addWorksheet(Enrichments, "Reactome C2")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP C0", x =(as.data.frame(EnrichR.FBs.C0$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF C0", x =(as.data.frame(EnrichR.FBs.C0$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC C0", x = (as.data.frame(EnrichR.FBs.C0$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome C0", x = (as.data.frame(EnrichR.FBs.C0$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP C2",  x = (as.data.frame(EnrichR.FBs.C2$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF C2",  x = (as.data.frame(EnrichR.FBs.C2$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC C2",  x = (as.data.frame(EnrichR.FBs.C2$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome C2", x = (as.data.frame(EnrichR.FBs.C2$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/FB C0 + C2 EnrichR.xlsx")

#Generating Enrichment bubble plots for marker genes in C0 (according to cut-off above)
EnrichR.FBs.C0.Top10<-as.data.frame(EnrichR.FBs.C0$GO_Biological_Process_2023)[1:10,]
frac <- EnrichR.FBs.C0.Top10$Overlap
EnrichR.FBs.C0.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(EnrichR.FBs.C0.Top10, aes(x=EnrichR.FBs.C0.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.C0.Top10$Adjusted.P.value, size=(EnrichR.FBs.C0.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/FB Cluster 0.BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

EnrichR.FBs.C0.CCTop10<-as.data.frame(EnrichR.FBs.C0$GO_Cellular_Component_2023)[1:10,]
frac <- EnrichR.FBs.C0.CCTop10$Overlap
EnrichR.FBs.C0.CCTop10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(EnrichR.FBs.C0.CCTop10, aes(x=EnrichR.FBs.C0.CCTop10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.C0.CCTop10$Adjusted.P.value, size=(EnrichR.FBs.C0.CCTop10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/FB Cluster 0.CC.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

EnrichR.FBs.C0.MF.Top10<-as.data.frame(EnrichR.FBs.C0$GO_Molecular_Function_2023)[1:10,]
frac <- EnrichR.FBs.C0.MF.Top10$Overlap
EnrichR.FBs.C0.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(EnrichR.FBs.C0.MF.Top10, aes(x=EnrichR.FBs.C0.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.C0.MF.Top10$Adjusted.P.value, size=(EnrichR.FBs.C0.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/FB Cluster 0.MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

EnrichR.FBs.C0.Reactome.Top10<-as.data.frame(EnrichR.FBs.C0$Reactome_2022)[1:10,]
frac <- EnrichR.FBs.C0.Reactome.Top10$Overlap
EnrichR.FBs.C0.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(EnrichR.FBs.C0.Reactome.Top10, aes(x=EnrichR.FBs.C0.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.C0.Reactome.Top10$Adjusted.P.value, size=(EnrichR.FBs.C0.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/FB Cluster 0.Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/FB Cluster 0.BP+CC+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, p4, labels = "AUTO", nrow = 1)
dev.off()

#Generating Enrichment bubble plots for marker genes in C2 (according to cut-off above)
EnrichR.FBs.C2.Top10<-as.data.frame(EnrichR.FBs.C0$GO_Biological_Process_2023)[1:10,]
frac <- EnrichR.FBs.C2.Top10$Overlap
EnrichR.FBs.C2.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(EnrichR.FBs.C2.Top10, aes(x=EnrichR.FBs.C2.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.C2.Top10$Adjusted.P.value, size=(EnrichR.FBs.C2.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/FB Cluster 2.BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

EnrichR.FBs.C2.CCTop10<-as.data.frame(EnrichR.FBs.C0$GO_Cellular_Component_2023)[1:10,]
frac <- EnrichR.FBs.C2.CCTop10$Overlap
EnrichR.FBs.C2.CCTop10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(EnrichR.FBs.C2.CCTop10, aes(x=EnrichR.FBs.C2.CCTop10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.C2.CCTop10$Adjusted.P.value, size=(EnrichR.FBs.C2.CCTop10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/FB Cluster 2.CC.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

EnrichR.FBs.C2.MF.Top10<-as.data.frame(EnrichR.FBs.C0$GO_Molecular_Function_2023)[1:10,]
frac <- EnrichR.FBs.C2.MF.Top10$Overlap
EnrichR.FBs.C2.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(EnrichR.FBs.C2.MF.Top10, aes(x=EnrichR.FBs.C2.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.C2.MF.Top10$Adjusted.P.value, size=(EnrichR.FBs.C2.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/FB Cluster 2.MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

EnrichR.FBs.C2.Reactome.Top10<-as.data.frame(EnrichR.FBs.C0$Reactome_2022)[1:10,]
frac <- EnrichR.FBs.C2.Reactome.Top10$Overlap
EnrichR.FBs.C2.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(EnrichR.FBs.C2.Reactome.Top10, aes(x=EnrichR.FBs.C2.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.C2.Reactome.Top10$Adjusted.P.value, size=(EnrichR.FBs.C2.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/FB Cluster 2.Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/FB Cluster 2.BP+CC+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, p4, labels = "AUTO", nrow = 1)
dev.off()


### Additional FB sub-cluster analysis - analysis for WT Fbs vs KO Fbs ####
Idents(ds_FBs) <- "ExpGroup"
FBs_markers.WTvsKO <- FindAllMarkers(ds_FBs, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_FB_WT vs KO.txt")

# Differential Gene Expression comparing WT and KO (non-pseudobulk): Stratified across time points
Idents(ds_FBs) <- "DPW"
table(ds_FBs$DPW)
ds_FBs.D0<-subset(ds_FBs, ident = "D0")
ds_FBs.D3<-subset(ds_FBs, ident = "D3")
ds_FBs.D7<-subset(ds_FBs, ident = "D7")

Idents(ds_FBs.D0) <- "ExpGroup"
FBs_markers.WTvsKO.D0 <- FindAllMarkers(ds_FBs.D0, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.D0, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.D0.txt")

Idents(ds_FBs.D3) <- "ExpGroup"
FBs_markers.WTvsKO.D3 <- FindAllMarkers(ds_FBs.D3, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.D3, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.D3.txt")

Idents(ds_FBs.D7) <- "ExpGroup"
FBs_markers.WTvsKO.D7 <- FindAllMarkers(ds_FBs.D7, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.D7, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.D7.txt")

###Making volcanoe plots for non-pseudobulk differential gene expression analysis over DPW
FBs_markers.WTvsKO.D0<-read.csv(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Nonpseudobulk Diff Exp Analysis_D0.csv")
FBs_markers.WTvsKO.D3<-read.csv(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Nonpseudobulk Diff Exp Analysis_D3.csv")
FBs_markers.WTvsKO.D7<-read.csv(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Nonpseudobulk Diff Exp Analysis_D7.csv")

nrow(FBs_markers.WTvsKO.D0) #3470 DEGs
nrow(FBs_markers.WTvsKO.D3) #1822 DEGs
nrow(FBs_markers.WTvsKO.D7) #721 DEGs

FBs_markers.WTvsKO.D0<- FBs_markers.WTvsKO.D0 %>%mutate(avg_log2FC = ifelse(cluster == "KO", avg_log2FC * -1, avg_log2FC))
FBs_markers.WTvsKO.D3<- FBs_markers.WTvsKO.D3 %>%mutate(avg_log2FC = ifelse(cluster == "KO", avg_log2FC * -1, avg_log2FC))
FBs_markers.WTvsKO.D7<- FBs_markers.WTvsKO.D7 %>%mutate(avg_log2FC = ifelse(cluster == "KO", avg_log2FC * -1, avg_log2FC))

FBs_markers.WTvsKO.D0.WT<-FBs_markers.WTvsKO.D0%>%filter(p_val_adj <  0.01 & avg_log2FC > 0.25 & pct.1 > 0.25)
FBs_markers.WTvsKO.D3.WT<-FBs_markers.WTvsKO.D3%>%filter(p_val_adj <  0.01 & avg_log2FC > 0.25 & pct.1 > 0.25)
FBs_markers.WTvsKO.D7.WT<-FBs_markers.WTvsKO.D7%>%filter(p_val_adj <  0.01 & avg_log2FC > 0.25 & pct.1 > 0.25)

FBs_markers.WTvsKO.D0.KO<-FBs_markers.WTvsKO.D0%>%filter(p_val_adj <  0.01 & avg_log2FC < -0.25 & pct.1 > 0.25)
FBs_markers.WTvsKO.D3.KO<-FBs_markers.WTvsKO.D3%>%filter(p_val_adj <  0.01 & avg_log2FC < -0.25 & pct.1 > 0.25)
FBs_markers.WTvsKO.D7.KO<-FBs_markers.WTvsKO.D7%>%filter(p_val_adj <  0.01 & avg_log2FC < -0.25 & pct.1 > 0.25)

FBs_markers.WTvsKO.D0<-rbind(FBs_markers.WTvsKO.D0.WT, FBs_markers.WTvsKO.D0.KO)
FBs_markers.WTvsKO.D3<-rbind(FBs_markers.WTvsKO.D3.WT, FBs_markers.WTvsKO.D3.KO)
FBs_markers.WTvsKO.D7<-rbind(FBs_markers.WTvsKO.D7.WT, FBs_markers.WTvsKO.D7.KO)

ggplot(data=FBs_markers.WTvsKO.D0, aes(x=avg_log2FC, y=-log10(p_val_adj), col=cluster)) + 
  theme(axis.title.x=element_text(size = 20), axis.title.y=element_text(size = 20), axis.text.x=element_text(size = 20), axis.text.y=element_text(size = 20), 
        legend.text = element_text(size=20), legend.title = element_blank(), panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = "grey", linewidth = 0.25),panel.grid.minor = element_line(color = "grey", linewidth = 0.25) )+
  geom_point(aes(colour=cluster, fill=cluster), size = 4, alpha = 1, shape = 21,colour = "black") + 
  geom_vline(xintercept=c(-0.25,0.25), col="red", size = 1.0) +   geom_hline(yintercept=-log10(0.01), col="red", size = 1.0)+
  xlab("Log2FC") + ylab("-log10(p.adj)")
nrow(FBs_markers.WTvsKO.D0.WT)
nrow(FBs_markers.WTvsKO.D0.KO)

ggplot(data=FBs_markers.WTvsKO.D3, aes(x=avg_log2FC, y=-log10(p_val_adj), col=cluster)) + 
  theme(axis.title.x=element_text(size = 20), axis.title.y=element_text(size = 20), axis.text.x=element_text(size = 20), axis.text.y=element_text(size = 20), 
        legend.text = element_text(size=20), legend.title = element_blank(), panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = "grey", linewidth = 0.25),panel.grid.minor = element_line(color = "grey", linewidth = 0.25) )+
  geom_point(aes(colour=cluster, fill=cluster), size = 4, alpha = 1, shape = 21,colour = "black") + 
  geom_vline(xintercept=c(-0.25,0.25), col="red", size = 1.0) +   geom_hline(yintercept=-log10(0.01), col="red", size = 1.0)+
  xlab("Log2FC") + ylab("-log10(p.adj)")
nrow(FBs_markers.WTvsKO.D3.WT)
nrow(FBs_markers.WTvsKO.D3.KO)

ggplot(data=FBs_markers.WTvsKO.D7, aes(x=avg_log2FC, y=-log10(p_val_adj), col=cluster)) + 
  theme(axis.title.x=element_text(size = 20), axis.title.y=element_text(size = 20), axis.text.x=element_text(size = 20), axis.text.y=element_text(size = 20), 
        legend.text = element_text(size=20), legend.title = element_blank(), panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = "grey", linewidth = 0.25),panel.grid.minor = element_line(color = "grey", linewidth = 0.25) )+
  geom_point(aes(colour=cluster, fill=cluster), size = 4, alpha = 1, shape = 21,colour = "black") + 
  geom_vline(xintercept=c(-0.25,0.25), col="red", size = 1.0) +   geom_hline(yintercept=-log10(0.01), col="red", size = 1.0)+
  xlab("Log2FC") + ylab("-log10(p.adj)")
nrow(FBs_markers.WTvsKO.D7.WT)
nrow(FBs_markers.WTvsKO.D7.KO)

###Enrichment of non-pseudobulk differential gene expression analysis across all FBs
FBs_markers.WTvsKO.D0<-read.csv(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Nonpseudobulk Diff Exp Analysis_D0.csv")
FBs_markers.WTvsKO.D3<-read.csv(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Nonpseudobulk Diff Exp Analysis_D3.csv")
FBs_markers.WTvsKO.D7<-read.csv(file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Nonpseudobulk Diff Exp Analysis_D7.csv")

FBs_markers.WTvsKO.D0<-FBs_markers.WTvsKO.D0%>%filter(p_val_adj <  0.01 & avg_log2FC > 0.25 & pct.1 > 0.25)
FBs_markers.WTvsKO.D3<-FBs_markers.WTvsKO.D3%>%filter(p_val_adj <  0.01 & avg_log2FC > 0.25 & pct.1 > 0.25)
FBs_markers.WTvsKO.D7<-FBs_markers.WTvsKO.D7%>%filter(p_val_adj <  0.01 & avg_log2FC > 0.25 & pct.1 > 0.25)

nrow(FBs_markers.WTvsKO) #357 DEGs
nrow(FBs_markers.WTvsKO.D0) #2175 DEGs
nrow(FBs_markers.WTvsKO.D3) #970 DEGs
nrow(FBs_markers.WTvsKO.D7) #147 DEGs

#non-stratified by DPW
FBs_markers.WTvsKO<-FBs_markers.WTvsKO%>%filter(p_val_adj <  0.01 & avg_log2FC > 0.25 & pct.1 > 0.25)
FBs_markers.WTvsKO.WT<-subset(FBs_markers.WTvsKO, cluster == "WT")
FBs_markers.WTvsKO.KO<-subset(FBs_markers.WTvsKO, cluster == "KO")

EnrichR.FBs.WTvsKO.WT<- enrichr(FBs_markers.WTvsKO.WT$gene, dbs)
EnrichR.FBs.WTvsKO.KO <- enrichr(FBs_markers.WTvsKO.KO$gene, dbs)

# non-stratified comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up WT")
addWorksheet(Enrichments, "MF Up WT")
addWorksheet(Enrichments, "CC Up WT")
addWorksheet(Enrichments, "Reactome Up WT")
addWorksheet(Enrichments, "BP Up KO")
addWorksheet(Enrichments, "MF Up KO")
addWorksheet(Enrichments, "CC Up KO")
addWorksheet(Enrichments, "Reactome Up KO")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up WT", x =(as.data.frame(EnrichR.FBs.WTvsKO.WT$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF Up WT", x =(as.data.frame(EnrichR.FBs.WTvsKO.WT$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC Up WT", x = (as.data.frame(EnrichR.FBs.WTvsKO.WT$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome Up WT", x = (as.data.frame(EnrichR.FBs.WTvsKO.WT$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.KO$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.KO$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.KO$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome Up KO", x = (as.data.frame(EnrichR.FBs.WTvsKO.KO$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/Not Stratified by DPW.xlsx")

#Generating Enrichment bubble plots for genes upregulated in WT
EnrichR.FBs.WTvsKO.WT.BP.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.WT$GO_Biological_Process_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.WT.BP.Top10$Overlap
EnrichR.FBs.WTvsKO.WT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(EnrichR.FBs.WTvsKO.WT.BP.Top10, aes(x=EnrichR.FBs.WTvsKO.WT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.WT.BP.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.WT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/No DPW split_WT vs KO.BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

EnrichR.FBs.WTvsKO.WT.CCTop10<-as.data.frame(EnrichR.FBs.WTvsKO.WT$GO_Cellular_Component_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.WT.CCTop10$Overlap
EnrichR.FBs.WTvsKO.WT.CCTop10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(EnrichR.FBs.WTvsKO.WT.CCTop10, aes(x=EnrichR.FBs.WTvsKO.WT.CCTop10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.WT.CCTop10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.WT.CCTop10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/No DPW split_WT vs KO.CC.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

EnrichR.FBs.WTvsKO.WT.MF.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.WT$GO_Molecular_Function_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.WT.MF.Top10$Overlap
EnrichR.FBs.WTvsKO.WT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(EnrichR.FBs.WTvsKO.WT.MF.Top10, aes(x=EnrichR.FBs.WTvsKO.WT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.WT.MF.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.WT.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/No DPW split_WT vs KO.MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

EnrichR.FBs.WTvsKO.WT.Reactome.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.WT$Reactome_2022)[1:10,]
frac <- EnrichR.FBs.WTvsKO.WT.Reactome.Top10$Overlap
EnrichR.FBs.WTvsKO.WT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(EnrichR.FBs.WTvsKO.WT.Reactome.Top10, aes(x=EnrichR.FBs.WTvsKO.WT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.WT.Reactome.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.WT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/No DPW split_WT vs KO.Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/No DPW split_WT vs KO.BP+CC+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, p4, labels = "AUTO", nrow = 1)
dev.off()

#Generating Enrichment bubble plots for genes upregulated in KO
EnrichR.FBs.WTvsKO.KO.BP.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.KO$GO_Biological_Process_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.KO.BP.Top10$Overlap
EnrichR.FBs.WTvsKO.KO.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(EnrichR.FBs.WTvsKO.KO.BP.Top10, aes(x=EnrichR.FBs.WTvsKO.KO.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.KO.BP.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.KO.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/No DPW split_WT vs KO.BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

EnrichR.FBs.WTvsKO.KO.CCTop10<-as.data.frame(EnrichR.FBs.WTvsKO.KO$GO_Cellular_Component_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.KO.CCTop10$Overlap
EnrichR.FBs.WTvsKO.KO.CCTop10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(EnrichR.FBs.WTvsKO.KO.CCTop10, aes(x=EnrichR.FBs.WTvsKO.KO.CCTop10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.KO.CCTop10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.KO.CCTop10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/No DPW split_WT vs KO.CC.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

EnrichR.FBs.WTvsKO.KO.MF.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.KO$GO_Molecular_Function_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.KO.MF.Top10$Overlap
EnrichR.FBs.WTvsKO.KO.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(EnrichR.FBs.WTvsKO.KO.MF.Top10, aes(x=EnrichR.FBs.WTvsKO.KO.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.KO.MF.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.KO.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/No DPW split_WT vs KO.MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

EnrichR.FBs.WTvsKO.KO.Reactome.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.KO$Reactome_2022)[1:10,]
frac <- EnrichR.FBs.WTvsKO.KO.Reactome.Top10$Overlap
EnrichR.FBs.WTvsKO.KO.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(EnrichR.FBs.WTvsKO.KO.Reactome.Top10, aes(x=EnrichR.FBs.WTvsKO.KO.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.KO.Reactome.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.KO.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/No DPW split_WT vs KO.Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/No DPW split_WT vs KO.BP+CC+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, p4, labels = "AUTO", nrow = 1)
dev.off()

#D0
FBs_markers.WTvsKO.D0<-FBs_markers.WTvsKO.D0%>%filter(p_val_adj <  0.01 & avg_log2FC > 0.25 & pct.1 > 0.25)
FBs_markers.WTvsKO.D0.WT<-subset(FBs_markers.WTvsKO.D0, cluster == "WT")
FBs_markers.WTvsKO.D0.KO<-subset(FBs_markers.WTvsKO.D0, cluster == "KO")

EnrichR.FBs.WTvsKO.D0.WT<- enrichr(FBs_markers.WTvsKO.D0.WT$gene, dbs)
EnrichR.FBs.WTvsKO.D0.KO <- enrichr(FBs_markers.WTvsKO.D0.KO$gene, dbs)

Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up WT")
addWorksheet(Enrichments, "MF Up WT")
addWorksheet(Enrichments, "CC Up WT")
addWorksheet(Enrichments, "Reactome Up WT")
addWorksheet(Enrichments, "BP Up KO")
addWorksheet(Enrichments, "MF Up KO")
addWorksheet(Enrichments, "CC Up KO")
addWorksheet(Enrichments, "Reactome Up KO")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up WT", x =(as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF Up WT", x =(as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC Up WT", x = (as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome Up WT", x = (as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.D0.KO$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.D0.KO$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.D0.KO$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome Up KO", x = (as.data.frame(EnrichR.FBs.WTvsKO.D0.KO$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.xlsx")

#Generating Enrichment bubble plots
EnrichR.FBs.WTvsKO.D0.WT.BP.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$GO_Biological_Process_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.WT.BP.Top10$Overlap
EnrichR.FBs.WTvsKO.D0.WT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(EnrichR.FBs.WTvsKO.D0.WT.BP.Top10, aes(x=EnrichR.FBs.WTvsKO.D0.WT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.WT.BP.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.WT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

EnrichR.FBs.WTvsKO.D0.WT.CCTop10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$GO_Cellular_Component_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.WT.CCTop10$Overlap
EnrichR.FBs.WTvsKO.D0.WT.CCTop10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(EnrichR.FBs.WTvsKO.D0.WT.CCTop10, aes(x=EnrichR.FBs.WTvsKO.D0.WT.CCTop10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.WT.CCTop10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.WT.CCTop10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.CC.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

EnrichR.FBs.WTvsKO.D0.WT.MF.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$GO_Molecular_Function_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.WT.MF.Top10$Overlap
EnrichR.FBs.WTvsKO.D0.WT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(EnrichR.FBs.WTvsKO.D0.WT.MF.Top10, aes(x=EnrichR.FBs.WTvsKO.D0.WT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.WT.MF.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.WT.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

EnrichR.FBs.WTvsKO.D0.WT.Reactome.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$Reactome_2022)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.WT.Reactome.Top10$Overlap
EnrichR.FBs.WTvsKO.D0.WT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(EnrichR.FBs.WTvsKO.D0.WT.Reactome.Top10, aes(x=EnrichR.FBs.WTvsKO.D0.WT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.WT.Reactome.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.WT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.BP+CC+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, p4, labels = "AUTO", nrow = 1)
dev.off()


#Generating Enrichment bubble plots for genes upregulated in KO
EnrichR.FBs.WTvsKO.D0.KO.BP.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$GO_Biological_Process_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.KO.BP.Top10$Overlap
EnrichR.FBs.WTvsKO.D0.KO.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(EnrichR.FBs.WTvsKO.D0.KO.BP.Top10, aes(x=EnrichR.FBs.WTvsKO.D0.KO.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.KO.BP.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.KO.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

EnrichR.FBs.WTvsKO.D0.KO.CCTop10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$GO_Cellular_Component_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.KO.CCTop10$Overlap
EnrichR.FBs.WTvsKO.D0.KO.CCTop10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(EnrichR.FBs.WTvsKO.D0.KO.CCTop10, aes(x=EnrichR.FBs.WTvsKO.D0.KO.CCTop10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.KO.CCTop10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.KO.CCTop10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.CC.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

EnrichR.FBs.WTvsKO.D0.KO.MF.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$GO_Molecular_Function_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.KO.MF.Top10$Overlap
EnrichR.FBs.WTvsKO.D0.KO.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(EnrichR.FBs.WTvsKO.D0.KO.MF.Top10, aes(x=EnrichR.FBs.WTvsKO.D0.KO.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.KO.MF.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.KO.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.WT$Reactome_2022)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10$Overlap
EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10, aes(x=EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.BP+CC+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, p4, labels = "AUTO", nrow = 1)
dev.off()


#Generating Enrichment bubble plots for genes upregulated in KO
EnrichR.FBs.WTvsKO.D0.KO.BP.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.KO$GO_Biological_Process_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.KO.BP.Top10$Overlap
EnrichR.FBs.WTvsKO.D0.KO.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(EnrichR.FBs.WTvsKO.D0.KO.BP.Top10, aes(x=EnrichR.FBs.WTvsKO.D0.KO.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.KO.BP.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.KO.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

EnrichR.FBs.WTvsKO.D0.KO.CCTop10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.KO$GO_Cellular_Component_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.KO.CCTop10$Overlap
EnrichR.FBs.WTvsKO.D0.KO.CCTop10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(EnrichR.FBs.WTvsKO.D0.KO.CCTop10, aes(x=EnrichR.FBs.WTvsKO.D0.KO.CCTop10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.KO.CCTop10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.KO.CCTop10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.CC.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

EnrichR.FBs.WTvsKO.D0.KO.MF.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.KO$GO_Molecular_Function_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.KO.MF.Top10$Overlap
EnrichR.FBs.WTvsKO.D0.KO.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(EnrichR.FBs.WTvsKO.D0.KO.MF.Top10, aes(x=EnrichR.FBs.WTvsKO.D0.KO.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.KO.MF.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.KO.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D0.KO$Reactome_2022)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10$Overlap
EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10, aes(x=EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D0.KO.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D0_WT vs KO.BP+CC+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, p4, labels = "AUTO", nrow = 1)
dev.off()

#D3
FBs_markers.WTvsKO.D3<-FBs_markers.WTvsKO.D3%>%filter(p_val_adj <  0.01 & avg_log2FC > 0.25 & pct.1 > 0.25)
FBs_markers.WTvsKO.D3.WT<-subset(FBs_markers.WTvsKO.D3.2, cluster == "WT")
FBs_markers.WTvsKO.D3.KO<-subset(FBs_markers.WTvsKO.D3.2, cluster == "KO")

EnrichR.FBs.WTvsKO.D3.WT<- enrichr(FBs_markers.WTvsKO.D3.WT$gene, dbs)
EnrichR.FBs.WTvsKO.D3.KO <- enrichr(FBs_markers.WTvsKO.D3.KO$gene, dbs)

Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up WT")
addWorksheet(Enrichments, "MF Up WT")
addWorksheet(Enrichments, "CC Up WT")
addWorksheet(Enrichments, "Reactome Up WT")
addWorksheet(Enrichments, "BP Up KO")
addWorksheet(Enrichments, "MF Up KO")
addWorksheet(Enrichments, "CC Up KO")
addWorksheet(Enrichments, "Reactome Up KO")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up WT", x =(as.data.frame(EnrichR.FBs.WTvsKO.D3.WT$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF Up WT", x =(as.data.frame(EnrichR.FBs.WTvsKO.D3.WT$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC Up WT", x = (as.data.frame(EnrichR.FBs.WTvsKO.D3.WT$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome Up WT", x = (as.data.frame(EnrichR.FBs.WTvsKO.D3.WT$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.D3.KO$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.D3.KO$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.D3.KO$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome Up KO", x = (as.data.frame(EnrichR.FBs.WTvsKO.D3.KO$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D3_WT vs KO.xlsx")

#Generating Enrichment bubble plots for genes upregulated in WT
EnrichR.FBs.WTvsKO.D3.WT.BP.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D3.WT$GO_Biological_Process_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D3.WT.BP.Top10$Overlap
EnrichR.FBs.WTvsKO.D3.WT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(EnrichR.FBs.WTvsKO.D3.WT.BP.Top10, aes(x=EnrichR.FBs.WTvsKO.D3.WT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D3.WT.BP.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D3.WT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D3_WT vs KO.BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

EnrichR.FBs.WTvsKO.D3.WT.CCTop10<-as.data.frame(EnrichR.FBs.WTvsKO.D3.WT$GO_Cellular_Component_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D3.WT.CCTop10$Overlap
EnrichR.FBs.WTvsKO.D3.WT.CCTop10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(EnrichR.FBs.WTvsKO.D3.WT.CCTop10, aes(x=EnrichR.FBs.WTvsKO.D3.WT.CCTop10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D3.WT.CCTop10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D3.WT.CCTop10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D3_WT vs KO.CC.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

EnrichR.FBs.WTvsKO.D3.WT.MF.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D3.WT$GO_Molecular_Function_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D3.WT.MF.Top10$Overlap
EnrichR.FBs.WTvsKO.D3.WT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(EnrichR.FBs.WTvsKO.D3.WT.MF.Top10, aes(x=EnrichR.FBs.WTvsKO.D3.WT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D3.WT.MF.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D3.WT.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D3_WT vs KO.MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

EnrichR.FBs.WTvsKO.D3.WT.Reactome.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D3.WT$Reactome_2022)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D3.WT.Reactome.Top10$Overlap
EnrichR.FBs.WTvsKO.D3.WT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(EnrichR.FBs.WTvsKO.D3.WT.Reactome.Top10, aes(x=EnrichR.FBs.WTvsKO.D3.WT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D3.WT.Reactome.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D3.WT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D3_WT vs KO.Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D3_WT vs KO.BP+CC+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, p4, labels = "AUTO", nrow = 1)
dev.off()

#Generating Enrichment bubble plots for genes upregulated in WT
EnrichR.FBs.WTvsKO.D3.KO.BP.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D3.KO$GO_Biological_Process_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D3.KO.BP.Top10$Overlap
EnrichR.FBs.WTvsKO.D3.KO.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(EnrichR.FBs.WTvsKO.D3.KO.BP.Top10, aes(x=EnrichR.FBs.WTvsKO.D3.KO.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D3.KO.BP.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D3.KO.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D3_WT vs KO.BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

EnrichR.FBs.WTvsKO.D3.KO.CCTop10<-as.data.frame(EnrichR.FBs.WTvsKO.D3.KO$GO_Cellular_Component_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D3.KO.CCTop10$Overlap
EnrichR.FBs.WTvsKO.D3.KO.CCTop10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(EnrichR.FBs.WTvsKO.D3.KO.CCTop10, aes(x=EnrichR.FBs.WTvsKO.D3.KO.CCTop10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D3.KO.CCTop10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D3.KO.CCTop10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D3_WT vs KO.CC.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

EnrichR.FBs.WTvsKO.D3.KO.MF.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D3.KO$GO_Molecular_Function_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D3.KO.MF.Top10$Overlap
EnrichR.FBs.WTvsKO.D3.KO.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(EnrichR.FBs.WTvsKO.D3.KO.MF.Top10, aes(x=EnrichR.FBs.WTvsKO.D3.KO.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D3.KO.MF.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D3.KO.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D3_WT vs KO.MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

EnrichR.FBs.WTvsKO.D3.KO.Reactome.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D3.KO$Reactome_2022)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D3.KO.Reactome.Top10$Overlap
EnrichR.FBs.WTvsKO.D3.KO.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(EnrichR.FBs.WTvsKO.D3.KO.Reactome.Top10, aes(x=EnrichR.FBs.WTvsKO.D3.KO.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D3.KO.Reactome.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D3.KO.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D3_WT vs KO.Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D3_WT vs KO.BP+CC+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, p4, labels = "AUTO", nrow = 1)
dev.off()

#D7
FBs_markers.WTvsKO.D7<-FBs_markers.WTvsKO.D7%>%filter(p_val_adj <  0.01 & avg_log2FC > 0.25 & pct.1 > 0.25)
FBs_markers.WTvsKO.D7.WT<-subset(FBs_markers.WTvsKO.D7, cluster == "WT")
FBs_markers.WTvsKO.D7.KO<-subset(FBs_markers.WTvsKO.D7, cluster == "KO")

EnrichR.FBs.WTvsKO.D7.WT<- enrichr(FBs_markers.WTvsKO.D7.WT$gene, dbs)
EnrichR.FBs.WTvsKO.D7.KO <- enrichr(FBs_markers.WTvsKO.D7.KO$gene, dbs)

Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up WT")
addWorksheet(Enrichments, "MF Up WT")
addWorksheet(Enrichments, "CC Up WT")
addWorksheet(Enrichments, "Reactome Up WT")
addWorksheet(Enrichments, "BP Up KO")
addWorksheet(Enrichments, "MF Up KO")
addWorksheet(Enrichments, "CC Up KO")
addWorksheet(Enrichments, "Reactome Up KO")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up WT", x =(as.data.frame(EnrichR.FBs.WTvsKO.D7.WT$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF Up WT", x =(as.data.frame(EnrichR.FBs.WTvsKO.D7.WT$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC Up WT", x = (as.data.frame(EnrichR.FBs.WTvsKO.D7.WT$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome Up WT", x = (as.data.frame(EnrichR.FBs.WTvsKO.D7.WT$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.D7.KO$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.D7.KO$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.D7.KO$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome Up KO", x = (as.data.frame(EnrichR.FBs.WTvsKO.D7.KO$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D7_WT vs KO.xlsx")

#Generating Enrichment bubble plots for genes upregulated in WT
EnrichR.FBs.WTvsKO.D7.WT.BP.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D7.WT$GO_Biological_Process_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D7.WT.BP.Top10$Overlap
EnrichR.FBs.WTvsKO.D7.WT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(EnrichR.FBs.WTvsKO.D7.WT.BP.Top10, aes(x=EnrichR.FBs.WTvsKO.D7.WT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D7.WT.BP.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D7.WT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D7_WT vs KO.BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

EnrichR.FBs.WTvsKO.D7.WT.CCTop10<-as.data.frame(EnrichR.FBs.WTvsKO.D7.WT$GO_Cellular_Component_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D7.WT.CCTop10$Overlap
EnrichR.FBs.WTvsKO.D7.WT.CCTop10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(EnrichR.FBs.WTvsKO.D7.WT.CCTop10, aes(x=EnrichR.FBs.WTvsKO.D7.WT.CCTop10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D7.WT.CCTop10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D7.WT.CCTop10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D7_WT vs KO.CC.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

#MF has no significance
EnrichR.FBs.WTvsKO.D7.WT.MF.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D7.WT$GO_Molecular_Function_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D7.WT.MF.Top10$Overlap
EnrichR.FBs.WTvsKO.D7.WT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(EnrichR.FBs.WTvsKO.D7.WT.MF.Top10, aes(x=EnrichR.FBs.WTvsKO.D7.WT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D7.WT.MF.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D7.WT.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D7_WT vs KO.MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

EnrichR.FBs.WTvsKO.D7.WT.Reactome.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D7.WT$Reactome_2022)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D7.WT.Reactome.Top10$Overlap
EnrichR.FBs.WTvsKO.D7.WT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(EnrichR.FBs.WTvsKO.D7.WT.Reactome.Top10, aes(x=EnrichR.FBs.WTvsKO.D7.WT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D7.WT.Reactome.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D7.WT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D7_WT vs KO.Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D7_WT vs KO.BP+CC+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p4, labels = "AUTO", nrow = 1)
dev.off()

#Generating Enrichment bubble plots for genes upregulated in KO
EnrichR.FBs.WTvsKO.D7.KO.BP.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D7.KO$GO_Biological_Process_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D7.KO.BP.Top10$Overlap
EnrichR.FBs.WTvsKO.D7.KO.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(EnrichR.FBs.WTvsKO.D7.KO.BP.Top10, aes(x=EnrichR.FBs.WTvsKO.D7.KO.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D7.KO.BP.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D7.KO.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D7_WT vs KO.BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

EnrichR.FBs.WTvsKO.D7.KO.CCTop10<-as.data.frame(EnrichR.FBs.WTvsKO.D7.KO$GO_Cellular_Component_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D7.KO.CCTop10$Overlap
EnrichR.FBs.WTvsKO.D7.KO.CCTop10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(EnrichR.FBs.WTvsKO.D7.KO.CCTop10, aes(x=EnrichR.FBs.WTvsKO.D7.KO.CCTop10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D7.KO.CCTop10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D7.KO.CCTop10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D7_WT vs KO.CC.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

#MF has no significance
EnrichR.FBs.WTvsKO.D7.KO.MF.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D7.KO$GO_Molecular_Function_2023)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D7.KO.MF.Top10$Overlap
EnrichR.FBs.WTvsKO.D7.KO.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(EnrichR.FBs.WTvsKO.D7.KO.MF.Top10, aes(x=EnrichR.FBs.WTvsKO.D7.KO.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D7.KO.MF.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D7.KO.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D7_WT vs KO.MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

EnrichR.FBs.WTvsKO.D7.KO.Reactome.Top10<-as.data.frame(EnrichR.FBs.WTvsKO.D7.KO$Reactome_2022)[1:10,]
frac <- EnrichR.FBs.WTvsKO.D7.KO.Reactome.Top10$Overlap
EnrichR.FBs.WTvsKO.D7.KO.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(EnrichR.FBs.WTvsKO.D7.KO.Reactome.Top10, aes(x=EnrichR.FBs.WTvsKO.D7.KO.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.FBs.WTvsKO.D7.KO.Reactome.Top10$Adjusted.P.value, size=(EnrichR.FBs.WTvsKO.D7.KO.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D7_WT vs KO.Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

png("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/D7_WT vs KO.BP+CC++MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, p4, labels = "AUTO", nrow = 1)
dev.off()


# Looking at the sub-cluster 0 and 2 for WT vs KO DEGs - not stratified by DPW
Idents(ds_FBs) <- "seurat_clusters"
table(ds_FBs$seurat_clusters)
ds_FBs.C0<-subset(ds_FBs, ident = "0")
ds_FBs.C2<-subset(ds_FBs, ident = "2")

Idents(ds_FBs.C0) <- "ExpGroup"
FBs_markers.WTvsKO.C0 <- FindAllMarkers(ds_FBs.C0, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.C0, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.Cluster 0.txt")

Idents(ds_FBs.C2) <- "ExpGroup"
FBs_markers.WTvsKO.C2 <- FindAllMarkers(ds_FBs.C2, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.C2, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.Cluster 2.txt")

#Enrichment of marker genes in Cluster 0 and Cluster 2: Upregulated in WT
FBs_markers.WTvsKO.C0<-FBs_markers.WTvsKO.C0%>%filter(p_val_adj <  0.01 & avg_log2FC > 0.25 & pct.1 > 0.25)

EnrichR.FBs.WTvsKO.WT<- enrichr(FBs_markers.WTvsKO.WT$gene, dbs)
EnrichR.FBs.WTvsKO.KO <- enrichr(FBs_markers.WTvsKO.KO$gene, dbs)

# non-stratified comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up WT")
addWorksheet(Enrichments, "MF Up WT")
addWorksheet(Enrichments, "CC Up WT")
addWorksheet(Enrichments, "Reactome Up WT")
addWorksheet(Enrichments, "BP Up KO")
addWorksheet(Enrichments, "MF Up KO")
addWorksheet(Enrichments, "CC Up KO")
addWorksheet(Enrichments, "Reactome Up KO")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up WT", x =(as.data.frame(EnrichR.FBs.WTvsKO.WT$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF Up WT", x =(as.data.frame(EnrichR.FBs.WTvsKO.WT$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC Up WT", x = (as.data.frame(EnrichR.FBs.WTvsKO.WT$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome Up WT", x = (as.data.frame(EnrichR.FBs.WTvsKO.WT$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.KO$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "MF Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.KO$GO_Molecular_Function_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "CC Up KO",  x = (as.data.frame(EnrichR.FBs.WTvsKO.KO$GO_Cellular_Component_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome Up KO", x = (as.data.frame(EnrichR.FBs.WTvsKO.KO$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk EnrichR/Not Stratified by DPW.xlsx")

### Additional FB sub-cluster analysis - analysis for WT Fbs vs KO Fbs: subcluster 0 and 2 ####
# Looking at the sub-cluster 0 and 2 for WT vs KO DEGs -  stratified by DPW
Idents(ds_FBs.C0) <- "DPW"
table(ds_FBs.C0$DPW)
ds_FBs.C0.D0<-subset(ds_FBs.C0, ident = "D0")
ds_FBs.C0.D3<-subset(ds_FBs.C0, ident = "D3")
ds_FBs.C0.D7<-subset(ds_FBs.C0, ident = "D7")

Idents(ds_FBs.C0.D0) <- "ExpGroup"
FBs_markers.WTvsKO.C0.D0 <- FindAllMarkers(ds_FBs.C0.D0, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.C0.D0, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.Cluster 0.D0.txt")

Idents(ds_FBs.C0.D3) <- "ExpGroup"
FBs_markers.WTvsKO.C0.D3 <- FindAllMarkers(ds_FBs.C0.D3, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.C0.D3, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.Cluster 0.D3.txt")

Idents(ds_FBs.C0.D7) <- "ExpGroup"
FBs_markers.WTvsKO.C0.D7 <- FindAllMarkers(ds_FBs.C0.D7, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.C0.D7, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.Cluster 0.D7.txt")

Idents(ds_FBs.C2) <- "DPW"
table(ds_FBs.C2$DPW)
ds_FBs.C2.D0<-subset(ds_FBs.C2, ident = "D0")
ds_FBs.C2.D3<-subset(ds_FBs.C2, ident = "D3")
ds_FBs.C2.D7<-subset(ds_FBs.C2, ident = "D7")

Idents(ds_FBs.C2.D0) <- "ExpGroup"
FBs_markers.WTvsKO.C2.D0 <- FindAllMarkers(ds_FBs.C2.D0, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.C2.D0, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.Cluster 2.D0.txt")

Idents(ds_FBs.C2.D3) <- "ExpGroup"
FBs_markers.WTvsKO.C2.D3 <- FindAllMarkers(ds_FBs.C2.D3, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.C2.D3, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.Cluster 2.D3.txt")

Idents(ds_FBs.C2.D7) <- "ExpGroup"
FBs_markers.WTvsKO.C2.D7 <- FindAllMarkers(ds_FBs.C2.D7, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25) 
write.csv(FBs_markers.WTvsKO.C2.D7, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/dim15res0.1_WTvsKO.Cluster 2.D7.txt")
# Additional Step: Enrichment of other FB clusters (1, 3-5) ####
# library(enrichR)
# listEnrichrSites()
# setEnrichrSite("Enrichr") # Human genes
# dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023","Reactome_2022")

fibroblast_markers.C1<-fibroblast_markers%>%filter(cluster == "1" & avg_log2FC >=1)
fibroblast_markers.C3<-fibroblast_markers%>%filter(cluster == "3" & avg_log2FC >=1)
fibroblast_markers.C4<-fibroblast_markers%>%filter(cluster == "4" & avg_log2FC >=1)
fibroblast_markers.C5<-fibroblast_markers%>%filter(cluster == "5" & avg_log2FC >=1)

# Cluster 1 enrichment
FB.C1.enriched <- enrichr(fibroblast_markers.C1$gene, dbs)

# Cluster 3 enrichment
FB.C3.enriched <- enrichr(fibroblast_markers.C3$gene, dbs)

# Cluster 4 enrichment
FB.C4.enriched <- enrichr(fibroblast_markers.C4$gene, dbs)

# Cluster 5 enrichment
FB.C5.enriched <- enrichr(fibroblast_markers.C5$gene, dbs)

#### Exporting the BP and Reactome Enrichment 
# Cluster 1, 3-5
Enrichment <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP C1")
addWorksheet(Enrichments, "CC C1")
addWorksheet(Enrichments, "MF C1")
addWorksheet(Enrichments, "Reactome C1")
addWorksheet(Enrichments, "BP C3")
addWorksheet(Enrichments, "CC C3")
addWorksheet(Enrichments, "MF C3")
addWorksheet(Enrichments, "Reactome C3")
addWorksheet(Enrichments, "BP C4")
addWorksheet(Enrichments, "CC C4")
addWorksheet(Enrichments, "MF C4")
addWorksheet(Enrichments, "Reactome C4")
addWorksheet(Enrichments, "BP C5")
addWorksheet(Enrichments, "CC C5")
addWorksheet(Enrichments, "MF C5")
addWorksheet(Enrichments, "Reactome C5")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP C1", x = as.data.frame(FB.C1.enriched$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC C1", x = as.data.frame(FB.C1.enriched$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF C1",  x = as.data.frame(FB.C1.enriched$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome C1", x = as.data.frame(FB.C1.enriched$Reactome_2022))
writeData(Enrichments, sheet = "BP C3", x = as.data.frame(FB.C3.enriched$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC C3", x = as.data.frame(FB.C3.enriched$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF C3",  x = as.data.frame(FB.C3.enriched$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome C3", x = as.data.frame(FB.C3.enriched$Reactome_2022))
writeData(Enrichments, sheet = "BP C4", x = as.data.frame(FB.C4.enriched$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC C4", x = as.data.frame(FB.C4.enriched$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF C4",  x = as.data.frame(FB.C4.enriched$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome C4", x = as.data.frame(FB.C4.enriched$Reactome_2022))
writeData(Enrichments, sheet = "BP C5", x = as.data.frame(FB.C5.enriched$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC C5", x = as.data.frame(FB.C5.enriched$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF C5",  x = as.data.frame(FB.C5.enriched$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome C5", x = as.data.frame(FB.C5.enriched$Reactome_2022))
# Export the file
saveWorkbook(Enrichments,  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/EnrichR_FB clusters 1 + 3-5.xlsx")

# Differential Expression of WT vs KO All FBs ####
table(dataset_fibroblast$Group) #WT vs KO
table(dataset_fibroblast$bulk.ident) #time-series analysis
table(dataset_fibroblast$orig.identity) #time-series analysis + Sex

# Getting "DEGs" (mnarker genes between KO and WT) by Running FindAllMarkers
#DS = dataset #NPB = non-pseudobulk
Idents(dataset_fibroblast) <- "Group"
DS_FBs.NPB.WTvsKO<- FindAllMarkers(dataset_fibroblast, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
write.csv(DS_FBs.NPB.WTvsKO, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk DEGs_WTvsKO FBs.txt") 
DS_FBs.NPB.WTvsKO<-read.xlsx("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk DEGs_WTvsKO FBs.xlsx") 

# cut-off 1: padj<0.01
DS_FBs.NPB.WTvsKO<-DS_FBs.NPB.WTvsKO%>%filter(p_val_adj < 0.01 )
DS_FBs.NPB.WTvsKO.WT<-DS_FBs.NPB.WTvsKO%>%filter(cluster == "WT" )
DS_FBs.NPB.WTvsKO.KO<-DS_FBs.NPB.WTvsKO%>%filter(cluster == "KO" )

#Enrichment via EnrichR
DS_FBs.NPB.WTvsKO.WT.enriched<-enrichr(DS_FBs.NPB.WTvsKO.WT$gene, dbs)
DS_FBs.NPB.WTvsKO.KO.enriched<-enrichr(DS_FBs.NPB.WTvsKO.KO$gene, dbs)

# cut-off 2: padj<0.0001
DS_FBs.NPB.WTvsKO<-DS_FBs.NPB.WTvsKO%>%filter(p_val_adj < 0.0001)
DS_FBs.NPB.WTvsKO.WT.C2<-DS_FBs.NPB.WTvsKO%>%filter(cluster == "WT")
DS_FBs.NPB.WTvsKO.KO.C2<-DS_FBs.NPB.WTvsKO%>%filter(cluster == "KO")

#Enrichment via EnrichR
DS_FBs.NPB.WTvsKO.WT.enriched.C2<-enrichr(DS_FBs.NPB.WTvsKO.WT.C2$gene, dbs)
DS_FBs.NPB.WTvsKO.KO.enriched.C2<-enrichr(DS_FBs.NPB.WTvsKO.KO.C2$gene, dbs)

# cut-off 3: padj<0.000001
DS_FBs.NPB.WTvsKO<-DS_FBs.NPB.WTvsKO%>%filter(p_val_adj < 0.000001)
DS_FBs.NPB.WTvsKO.WT.C3<-DS_FBs.NPB.WTvsKO%>%filter(cluster == "WT")
DS_FBs.NPB.WTvsKO.KO.C3<-DS_FBs.NPB.WTvsKO%>%filter(cluster == "KO")

#Enrichment via EnrichR
DS_FBs.NPB.WTvsKO.WT.enriched.C3<-enrichr(DS_FBs.NPB.WTvsKO.WT.C3$gene, dbs)
DS_FBs.NPB.WTvsKO.KO.enriched.C3<-enrichr(DS_FBs.NPB.WTvsKO.KO.C3$gene, dbs)

#export: cutoff 1 p < 0.01
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP WT")
addWorksheet(Enrichments, "CC WT")
addWorksheet(Enrichments, "MF WT")
addWorksheet(Enrichments, "Reactome WT")
addWorksheet(Enrichments, "BP KO")
addWorksheet(Enrichments, "CC KO")
addWorksheet(Enrichments, "MF KO")
addWorksheet(Enrichments, "Reactome KO")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP WT", x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC WT", x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF WT",  x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome WT", x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched$Reactome_2022))
writeData(Enrichments, sheet = "BP KO", x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC KO", x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF KO",  x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome KO", x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched$Reactome_2022))
# Export the file
saveWorkbook(Enrichments,  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/EnrichR_Non-pseudobulk_WT vs KO.C1.xlsx")

#export: cutoff 2 p < 0.0001
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP WT")
addWorksheet(Enrichments, "CC WT")
addWorksheet(Enrichments, "MF WT")
addWorksheet(Enrichments, "Reactome WT")
addWorksheet(Enrichments, "BP KO")
addWorksheet(Enrichments, "CC KO")
addWorksheet(Enrichments, "MF KO")
addWorksheet(Enrichments, "Reactome KO")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP WT", x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched.C2$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC WT", x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched.C2$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF WT",  x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched.C2$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome WT", x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched.C2$Reactome_2022))
writeData(Enrichments, sheet = "BP KO", x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched.C2$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC KO", x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched.C2$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF KO",  x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched.C2$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome KO", x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched.C2$Reactome_2022))
# Export the file
saveWorkbook(Enrichments,  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/EnrichR_Non-pseudobulk_WT vs KO.C2.xlsx")

#export: cutoff 3 p < 0.000001
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP WT")
addWorksheet(Enrichments, "CC WT")
addWorksheet(Enrichments, "MF WT")
addWorksheet(Enrichments, "Reactome WT")
addWorksheet(Enrichments, "BP KO")
addWorksheet(Enrichments, "CC KO")
addWorksheet(Enrichments, "MF KO")
addWorksheet(Enrichments, "Reactome KO")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP WT", x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched.C3$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC WT", x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched.C3$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF WT",  x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched.C3$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome WT", x = as.data.frame(DS_FBs.NPB.WTvsKO.WT.enriched.C3$Reactome_2022))
writeData(Enrichments, sheet = "BP KO", x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched.C3$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC KO", x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched.C3$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF KO",  x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched.C3$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome KO", x = as.data.frame(DS_FBs.NPB.WTvsKO.KO.enriched.C3$Reactome_2022))
# Export the file
saveWorkbook(Enrichments,  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/EnrichR_Non-pseudobulk_WT vs KO.C3.xlsx")

# Differential Expression of FB cluster 0 (MFAP5+) between WT and KO - not stratified by DPW ####
table(dataset_fibroblast.C0$Group) #comparing WT and KO
table(dataset_fibroblast$bulk.ident) #time-series analysis

# Getting "DEGs" (mnarker genes between KO and WT) by Running FindAllMarkers
#DS = dataset #NPB = non-pseudobulk
Idents(dataset_fibroblast.C0) <- "Group"
DS_FBs.C0.NPB.WTvsKO<- FindAllMarkers(dataset_fibroblast.C0, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
write.csv(DS_FBs.C0.NPB.WTvsKO, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk DEGs_WTvsKO C0 FBs.txt") 
DS_FBs.C0.NPB.WTvsKO<-read.xlsx("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk DEGs_WTvsKO C0 FBs.xlsx") 

# cut-off 1: padj<0.01
DS_FBs.C0.NPB.WTvsKO<-DS_FBs.C0.NPB.WTvsKO%>%filter(p_val_adj < 0.01 )
DS_FBs.C0.NPB.WTvsKO.WT<-DS_FBs.C0.NPB.WTvsKO%>%filter(cluster == "WT" )
DS_FBs.C0.NPB.WTvsKO.KO<-DS_FBs.C0.NPB.WTvsKO%>%filter(cluster == "KO" )

#Enrichment via EnrichR
DS_FBs.CO.NPB.WTvsKO.WT.enriched.C1<-enrichr(DS_FBs.C0.NPB.WTvsKO.WT$gene, dbs)
DS_FBs.CO.NPB.WTvsKO.KO.enriched.C1<-enrichr(DS_FBs.C0.NPB.WTvsKO.KO$gene, dbs)

# cut-off 2: padj<0.0001
DS_FBs.C0.NPB.WTvsKO.C2<-DS_FBs.C0.NPB.WTvsKO%>%filter(p_val_adj < 0.0001 )
DS_FBs.C0.NPB.WTvsKO.WT.C2<-DS_FBs.C0.NPB.WTvsKO.C2%>%filter(cluster == "WT" )
DS_FBs.C0.NPB.WTvsKO.KO.C2<-DS_FBs.C0.NPB.WTvsKO.C2%>%filter(cluster == "KO" )

#Enrichment via EnrichR
DS_FBs.CO.NPB.WTvsKO.WT.enriched.C2<-enrichr(DS_FBs.C0.NPB.WTvsKO.WT.C2$gene, dbs)
DS_FBs.CO.NPB.WTvsKO.KO.enriched.C2<-enrichr(DS_FBs.C0.NPB.WTvsKO.KO.C2$gene, dbs)

# cut-off 3: padj<0.000001
DS_FBs.C0.NPB.WTvsKO.C3<-DS_FBs.C0.NPB.WTvsKO%>%filter(p_val_adj < 0.000001 )
DS_FBs.C0.NPB.WTvsKO.WT.C3<-DS_FBs.C0.NPB.WTvsKO.C3%>%filter(cluster == "WT" )
DS_FBs.C0.NPB.WTvsKO.KO.C3<-DS_FBs.C0.NPB.WTvsKO.C3%>%filter(cluster == "KO" )

#Enrichment via EnrichR
DS_FBs.CO.NPB.WTvsKO.WT.enriched.C3<-enrichr(DS_FBs.C0.NPB.WTvsKO.WT.C3$gene, dbs)
DS_FBs.CO.NPB.WTvsKO.KO.enriched.C3<-enrichr(DS_FBs.C0.NPB.WTvsKO.KO.C3$gene, dbs)

#export: cutoff 1 p < 0.01
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP WT")
addWorksheet(Enrichments, "CC WT")
addWorksheet(Enrichments, "MF WT")
addWorksheet(Enrichments, "Reactome WT")
addWorksheet(Enrichments, "BP KO")
addWorksheet(Enrichments, "CC KO")
addWorksheet(Enrichments, "MF KO")
addWorksheet(Enrichments, "Reactome KO")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP WT", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C1$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC WT", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C1$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF WT",  x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C1$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome WT", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C1$Reactome_2022))
writeData(Enrichments, sheet = "BP KO", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C1$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC KO", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C1$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF KO",  x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C1$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome KO", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C1$Reactome_2022))
# Export the file
saveWorkbook(Enrichments,  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/EnrichR_Non-pseudobulk_WT vs KO.FB C0.C1.xlsx")

#export: cutoff 2 p < 0.0001
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP WT")
addWorksheet(Enrichments, "CC WT")
addWorksheet(Enrichments, "MF WT")
addWorksheet(Enrichments, "Reactome WT")
addWorksheet(Enrichments, "BP KO")
addWorksheet(Enrichments, "CC KO")
addWorksheet(Enrichments, "MF KO")
addWorksheet(Enrichments, "Reactome KO")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP WT", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C2$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC WT", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C2$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF WT",  x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C2$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome WT", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C2$Reactome_2022))
writeData(Enrichments, sheet = "BP KO", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C2$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC KO", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C2$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF KO",  x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C2$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome KO", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C2$Reactome_2022))
# Export the file
saveWorkbook(Enrichments,  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/EnrichR_Non-pseudobulk_WT vs KO.FB C0.C2.xlsx")

#export: cutoff 3 p < 0.000001
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP WT")
addWorksheet(Enrichments, "CC WT")
addWorksheet(Enrichments, "MF WT")
addWorksheet(Enrichments, "Reactome WT")
addWorksheet(Enrichments, "BP KO")
addWorksheet(Enrichments, "CC KO")
addWorksheet(Enrichments, "MF KO")
addWorksheet(Enrichments, "Reactome KO")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP WT", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C3$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC WT", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C3$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF WT",  x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C3$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome WT", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.WT.enriched.C3$Reactome_2022))
writeData(Enrichments, sheet = "BP KO", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C3$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC KO", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C3$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF KO",  x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C3$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome KO", x = as.data.frame(DS_FBs.CO.NPB.WTvsKO.KO.enriched.C3$Reactome_2022))
# Export the file
saveWorkbook(Enrichments,  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/EnrichR_Non-pseudobulk_WT vs KO.FB C0.C3.xlsx")


# Differential Expression of FB cluster 0 (MFAP5+) between WT and KO - stratified by DPW ####
table(dataset_fibroblast$bulk.ident) #time-series analysis

# Getting "DEGs" (mnarker genes between KO and WT) by Running FindAllMarkers
#DS = dataset #NPB = non-pseudobulk
Idents(dataset_fibroblast.C0) <- "bulk.ident"
DS_FBs.C0.NPB.WTvsKO.DPW<- FindAllMarkers(dataset_fibroblast.C0, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
write.csv(DS_FBs.C0.NPB.WTvsKO.DPW, file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk DEGs_WTvsKO C0 FBs.DPW.txt") 
DS_FBs.C0.NPB.WTvsKO.DPW<-read.xlsx("C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/Non-pseudobulk DEGs_WTvsKO C0 FBs.DPW.xlsx") 

# cut-off 1: padj<0.01
DS_FBs.C0.NPB.WTvsKO.DPW<-DS_FBs.C0.NPB.WTvsKO.DPW%>%filter(p_val_adj < 0.01 )
DS_FBs.C0.NPB.WTvsKO.DPW.KOD0<-DS_FBs.C0.NPB.WTvsKO.DPW%>%filter(cluster == "KOD0" )
DS_FBs.C0.NPB.WTvsKO.DPW.KOD3<-DS_FBs.C0.NPB.WTvsKO.DPW%>%filter(cluster == "KOD3" )
DS_FBs.C0.NPB.WTvsKO.DPW.KOD7<-DS_FBs.C0.NPB.WTvsKO.DPW%>%filter(cluster == "KOD7" )
DS_FBs.C0.NPB.WTvsKO.DPW.WTD0<-DS_FBs.C0.NPB.WTvsKO.DPW%>%filter(cluster == "WTD3" )
DS_FBs.C0.NPB.WTvsKO.DPW.WTD3<-DS_FBs.C0.NPB.WTvsKO.DPW%>%filter(cluster == "WTD0" )
DS_FBs.C0.NPB.WTvsKO.DPW.WTD7<-DS_FBs.C0.NPB.WTvsKO.DPW%>%filter(cluster == "WTD3" )

#Enrichment via EnrichR
DS_FBs.C0.NPB.WTvsKO.DPW.KOD0.enriched.C1<-enrichr(DS_FBs.C0.NPB.WTvsKO.DPW.KOD0$gene, dbs)
DS_FBs.C0.NPB.WTvsKO.DPW.KOD3.enriched.C1<-enrichr(DS_FBs.C0.NPB.WTvsKO.DPW.KOD3$gene, dbs)
DS_FBs.C0.NPB.WTvsKO.DPW.KOD7.enriched.C1<-enrichr(DS_FBs.C0.NPB.WTvsKO.DPW.KOD7$gene, dbs)
DS_FBs.C0.NPB.WTvsKO.DPW.WTD0.enriched.C1<-enrichr(DS_FBs.C0.NPB.WTvsKO.DPW.WTD0$gene, dbs)
DS_FBs.C0.NPB.WTvsKO.DPW.WTD3.enriched.C1<-enrichr(DS_FBs.C0.NPB.WTvsKO.DPW.WTD3$gene, dbs)
DS_FBs.C0.NPB.WTvsKO.DPW.WTD7.enriched.C1<-enrichr(DS_FBs.C0.NPB.WTvsKO.DPW.WTD7$gene, dbs)

#export: cutoff 1 p < 0.01
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP WTD0")
addWorksheet(Enrichments, "CC WTD0")
addWorksheet(Enrichments, "MF WTD0")
addWorksheet(Enrichments, "Reactome WTD0")
addWorksheet(Enrichments, "BP WTD3")
addWorksheet(Enrichments, "CC WTD3")
addWorksheet(Enrichments, "MF WTD3")
addWorksheet(Enrichments, "Reactome WTD3")
addWorksheet(Enrichments, "BP WTD7")
addWorksheet(Enrichments, "CC WTD7")
addWorksheet(Enrichments, "MF WTD7")
addWorksheet(Enrichments, "Reactome WTD7")
addWorksheet(Enrichments, "BP KOD0")
addWorksheet(Enrichments, "CC KOD0")
addWorksheet(Enrichments, "MF KOD0")
addWorksheet(Enrichments, "Reactome KOD0")
addWorksheet(Enrichments, "BP KOD3")
addWorksheet(Enrichments, "CC KOD3")
addWorksheet(Enrichments, "MF KOD3")
addWorksheet(Enrichments, "Reactome KOD3")
addWorksheet(Enrichments, "BP KOD7")
addWorksheet(Enrichments, "CC KOD7")
addWorksheet(Enrichments, "MF KOD7")
addWorksheet(Enrichments, "Reactome KOD7")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP WTD0", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD0.enriched.C1$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC WTD0", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD0.enriched.C1$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF WTD0",  x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD0.enriched.C1$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome WTD0", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD0.enriched.C1$Reactome_2022))
writeData(Enrichments, sheet = "BP WTD3", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD3.enriched.C1$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC WTD3", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD3.enriched.C1$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF WTD3",  x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD3.enriched.C1$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome WTD3", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD3.enriched.C1$Reactome_2022))
writeData(Enrichments, sheet = "BP WTD7", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD7.enriched.C1$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC WTD7", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD7.enriched.C1$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF WTD7",  x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD7.enriched.C1$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome WTD7", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.WTD7.enriched.C1$Reactome_2022))
writeData(Enrichments, sheet = "BP KOD0", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD0.enriched.C1$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC KOD0", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD0.enriched.C1$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF KOD0",  x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD0.enriched.C1$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome KOD0", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD0.enriched.C1$Reactome_2022))
writeData(Enrichments, sheet = "BP KOD3", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD3.enriched.C1$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC KOD3", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD3.enriched.C1$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF KOD3",  x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD3.enriched.C1$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome KOD3", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD3.enriched.C1$Reactome_2022))
writeData(Enrichments, sheet = "BP KOD7", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD7.enriched.C1$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "CC KOD7", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD7.enriched.C1$GO_Cellular_Component_2023))
writeData(Enrichments, sheet = "MF KOD7",  x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD7.enriched.C1$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome KOD7", x = as.data.frame(DS_FBs.C0.NPB.WTvsKO.DPW.KOD7.enriched.C1$Reactome_2022))

# Export the file
saveWorkbook(Enrichments,  file = "C:/Users/chenh/Desktop/sc-RNA Seq analysis/Merged Analysis/WT+KO Analysis/FB subpop/EnrichR_Non-pseudobulk_WT vs KO.FB C0 time-analysis.C1.xlsx")



