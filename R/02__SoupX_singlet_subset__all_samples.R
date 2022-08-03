library(dplyr)
library(Seurat)
library(SoupX)
library(data.table)
library(cowplot)
library(ggplot2)
set.seed(123)



#===============================
#============ SoupX ============
#===============================
# Load 10x datasets and create Seurat objects
RK1_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK1_4_comb/RK1_4_comb/outs/filtered_feature_bc_matrix")
RK3_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK3_4_comb/outs/filtered_feature_bc_matrix")
RK5_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK5_4_comb/outs/filtered_feature_bc_matrix")
RK7_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK7_4_comb/outs/filtered_feature_bc_matrix")
RK9_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK9_4_comb/outs/filtered_feature_bc_matrix")
RK16_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK16_4_comb/outs/filtered_feature_bc_matrix")
RK43_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK43_4_comb/RK43_4_comb/outs/filtered_feature_bc_matrix")
RK45_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK45_4_comb/RK45_4_comb/outs/filtered_feature_bc_matrix")
RK46_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK46_4_comb/RK46_4_comb/outs/filtered_feature_bc_matrix")
RK57_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK57_4_comb/RK57_4_comb/outs/filtered_feature_bc_matrix")
RK58_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK58_4_comb/RK58_4_comb/outs/filtered_feature_bc_matrix")
RK61_4.data <- Read10X(data.dir = "/~/cellranger_count/object_comb/RK61_4_comb/RK61_4_comb/outs/filtered_feature_bc_matrix")

RK1_4 <- CreateSeuratObject(counts = RK1_4.data, project = "RK1_4")
RK3_4 <- CreateSeuratObject(counts = RK3_4.data, project = "RK3_4")
RK5_4 <- CreateSeuratObject(counts = RK5_4.data, project = "RK5_4")
RK7_4 <- CreateSeuratObject(counts = RK7_4.data, project = "RK7_4")
RK9_4 <- CreateSeuratObject(counts = RK9_4.data, project = "RK9_4")
RK16_4 <- CreateSeuratObject(counts = RK16_4.data, project = "RK16_4")
RK43_4 <- CreateSeuratObject(counts = RK43_4.data, project = "RK43_4")
RK45_4 <- CreateSeuratObject(counts = RK45_4.data, project = "RK45_4")
RK46_4 <- CreateSeuratObject(counts = RK46_4.data, project = "RK46_4")
RK57_4 <- CreateSeuratObject(counts = RK57_4.data, project = "RK57_4")
RK58_4 <- CreateSeuratObject(counts = RK58_4.data, project = "RK58_4")
RK61_4 <- CreateSeuratObject(counts = RK61_4.data, project = "RK61_4")

# Load 10X data and estimate soup profile
sc_RK1_4 = load10X("/~/cellranger_count/object_comb/RK1_4_comb/RK1_4_comb/outs/") 
sc_RK3_4 = load10X("/~/cellranger_count/object_comb/RK3_4_comb/outs/") 
sc_RK5_4 = load10X("/~/cellranger_count/object_comb/RK5_4_comb/outs/") 
sc_RK7_4 = load10X("/~/cellranger_count/object_comb/RK7_4_comb/outs/") 
sc_RK9_4 = load10X("/~/cellranger_count/object_comb/RK9_4_comb/outs/") 
sc_RK16_4 = load10X("/~/cellranger_count/object_comb/RK16_4_comb/outs/") 
sc_RK43_4 = load10X("/~/cellranger_count/object_comb/RK43_4_comb/RK43_4_comb/outs/") 
sc_RK45_4 = load10X("/~/cellranger_count/object_comb/RK45_4_comb/RK45_4_comb/outs/") 
sc_RK46_4 = load10X("/~/cellranger_count/object_comb/RK46_4_comb/RK46_4_comb/outs/") 
sc_RK57_4 = load10X("/~/cellranger_count/object_comb/RK57_4_comb/RK57_4_comb/outs/") 
sc_RK58_4 = load10X("/~/cellranger_count/object_comb/RK58_4_comb/RK58_4_comb/outs/") 
sc_RK61_4 = load10X("/~/cellranger_count/object_comb/RK61_4_comb/RK61_4_comb/outs/") 

# Estimate rho
sc_RK1_4 = autoEstCont(sc_RK1_4) #Estimated global rho of 0.21
sc_RK3_4 = autoEstCont(sc_RK3_4) #Estimated global rho of 0.18
sc_RK5_4 = autoEstCont(sc_RK5_4) #Estimated global rho of 0.31
sc_RK7_4 = autoEstCont(sc_RK7_4) #Estimated global rho of 0.36
sc_RK9_4 = autoEstCont(sc_RK9_4) #Estimated global rho of 0.38
sc_RK16_4 = autoEstCont(sc_RK16_4) #Estimated global rho of 0.14
sc_RK43_4 = autoEstCont(sc_RK43_4) #Estimated global rho of 0.23
sc_RK45_4 = autoEstCont(sc_RK45_4) #Estimated global rho of 0.38
sc_RK46_4 = autoEstCont(sc_RK46_4) #Estimated global rho of 0.20
sc_RK57_4 = autoEstCont(sc_RK57_4) #Estimated global rho of 0.22
sc_RK58_4 = autoEstCont(sc_RK58_4) #Estimated global rho of 0.32
sc_RK61_4 = autoEstCont(sc_RK61_4) #Estimated global rho of 0.26

# Clean the data
out_RK1_4 = adjustCounts(sc_RK1_4, roundToInt=T)
out_RK3_4 = adjustCounts(sc_RK3_4, roundToInt=T)
out_RK5_4 = adjustCounts(sc_RK5_4, roundToInt=T)
out_RK7_4 = adjustCounts(sc_RK7_4, roundToInt=T)
out_RK9_4 = adjustCounts(sc_RK9_4, roundToInt=T)
out_RK16_4 = adjustCounts(sc_RK16_4, roundToInt=T)
out_RK43_4 = adjustCounts(sc_RK43_4, roundToInt=T)
out_RK45_4 = adjustCounts(sc_RK45_4, roundToInt=T)
out_RK46_4 = adjustCounts(sc_RK46_4, roundToInt=T)
out_RK57_4 = adjustCounts(sc_RK57_4, roundToInt=T)
out_RK58_4 = adjustCounts(sc_RK58_4, roundToInt=T)
out_RK61_4 = adjustCounts(sc_RK61_4, roundToInt=T)

# Overwrite Seurat object count matrix with ambient RNA-cleaned, inter-rounded matrix from SoupX
RK1_4@assays$RNA@counts@x <- out_RK1_4@x
RK3_4@assays$RNA@counts@x <- out_RK3_4@x
RK5_4@assays$RNA@counts@x <- out_RK5_4@x
RK7_4@assays$RNA@counts@x <- out_RK7_4@x
RK9_4@assays$RNA@counts@x <- out_RK9_4@x
RK16_4@assays$RNA@counts@x <- out_RK16_4@x
RK43_4@assays$RNA@counts@x <- out_RK43_4@x
RK45_4@assays$RNA@counts@x <- out_RK45_4@x
RK46_4@assays$RNA@counts@x <- out_RK46_4@x
RK57_4@assays$RNA@counts@x <- out_RK57_4@x
RK58_4@assays$RNA@counts@x <- out_RK58_4@x
RK61_4@assays$RNA@counts@x <- out_RK61_4@x

# Overwrite Seurat object data matrix with ambient RNA-cleaned, inter-rounded matrix from SoupX
RK1_4@assays$RNA@data@x <- out_RK1_4@x
RK3_4@assays$RNA@data@x <- out_RK3_4@x
RK5_4@assays$RNA@data@x <- out_RK5_4@x
RK7_4@assays$RNA@data@x <- out_RK7_4@x
RK9_4@assays$RNA@data@x <- out_RK9_4@x
RK16_4@assays$RNA@data@x <- out_RK16_4@x
RK43_4@assays$RNA@data@x <- out_RK43_4@x
RK45_4@assays$RNA@data@x <- out_RK45_4@x
RK46_4@assays$RNA@data@x <- out_RK46_4@x
RK57_4@assays$RNA@data@x <- out_RK57_4@x
RK58_4@assays$RNA@data@x <- out_RK58_4@x
RK61_4@assays$RNA@data@x <- out_RK61_4@x

# Remove clutter
rm(RK1_4.data,RK3_4.data,RK5_4.data,RK7_4.data,RK9_4.data,RK16_4.data,RK43_4.data,RK45_4.data,RK46_4.data,RK57_4.data,RK58_4.data,RK61_4.data)
rm(sc_RK1_4,sc_RK3_4,sc_RK5_4,sc_RK7_4,sc_RK9_4,sc_RK16_4,sc_RK43_4,sc_RK45_4,sc_RK46_4,sc_RK57_4,sc_RK58_4,sc_RK61_4)
rm(out_RK1_4,out_RK3_4,out_RK5_4,out_RK7_4,out_RK9_4,out_RK16_4,out_RK43_4,out_RK45_4,out_RK46_4,out_RK57_4,out_RK58_4,out_RK61_4)

# Merge Seurat objects
object <- merge(RK1_4, y = c(RK3_4, RK5_4, RK7_4, RK9_4, RK16_4, RK43_4, RK45_4, RK46_4, RK57_4, RK58_4, RK61_4), 
                     add.cell.ids = c("RK1_4", "RK3_4", "RK5_4", "RK7_4", "RK9_4", "RK16_4", "RK43_4", "RK45_4", "RK46_4", "RK57_4", "RK58_4", "RK61_4"),
                     project = "object_comb")
object #25399 features across 269668 samples within 1 assay



#=========================================
#============ SUBSET SINGLETS ============
#=========================================
# Import singlet barcodes for all the samples (DoubletFinder was run on every individual sample)
singlets_RK1_4 <- readRDS('/~/seurat/RK1_4_comb_soupX/singlets.rds')
singlets_RK3_4 <- readRDS('/~/seurat/RK3_4_comb_soupX/singlets.rds')
singlets_RK5_4 <- readRDS('/~/seurat/RK5_4_comb_soupX/singlets.rds')
singlets_RK7_4 <- readRDS('/~/seurat/RK7_4_comb_soupX/singlets.rds')
#singlets_RK9_4 <- readRDS('/~/seurat/RK9_4_comb_soupX/singlets.rds')
allcells_RK9_4 <- readRDS('/~/seurat/RK9_4_comb_soupX/barcodes.rds')
singlets_RK16_4 <- readRDS('/~/seurat/RK16_4_comb_soupX/singlets.rds')
singlets_RK43_4 <- readRDS('/~/seurat/RK43_4_comb_soupX/singlets.rds')
singlets_RK45_4 <- readRDS('/~/seurat/RK45_4_comb_soupX/singlets.rds')
singlets_RK46_4 <- readRDS('/~/seurat/RK46_4_comb_soupX/singlets.rds')
singlets_RK57_4 <- readRDS('/~/seurat/RK57_4_comb_soupX/singlets.rds')
singlets_RK58_4 <- readRDS('/~/seurat/RK58_4_comb_soupX/singlets.rds')
singlets_RK61_4 <- readRDS('/~/seurat/RK61_4_comb_soupX/singlets.rds')

# Add sample name as prefix to barcodes
singlets_RK1_4 <- paste0("RK1_4_",singlets_RK1_4)
singlets_RK3_4 <- paste0("RK3_4_",singlets_RK3_4)
singlets_RK5_4 <- paste0("RK5_4_",singlets_RK5_4)
singlets_RK7_4 <- paste0("RK7_4_",singlets_RK7_4)
#singlets_RK9_4 <- paste0("RK9_4_",singlets_RK9_4)
allcells_RK9_4 <- paste0("RK9_4_",allcells_RK9_4)
singlets_RK16_4 <- paste0("RK16_4_",singlets_RK16_4)
singlets_RK43_4 <- paste0("RK43_4_",singlets_RK43_4)
singlets_RK45_4 <- paste0("RK45_4_",singlets_RK45_4)
singlets_RK46_4 <- paste0("RK46_4_",singlets_RK46_4)
singlets_RK57_4 <- paste0("RK57_4_",singlets_RK57_4)
singlets_RK58_4 <- paste0("RK58_4_",singlets_RK58_4)
singlets_RK61_4 <- paste0("RK61_4_",singlets_RK61_4)

# Subset object
cells.to.subset <- c(singlets_RK1_4, singlets_RK3_4, singlets_RK5_4, singlets_RK7_4, #singlets_RK9_4, 
                     allcells_RK9_4, 
                     singlets_RK16_4, singlets_RK43_4, singlets_RK45_4, singlets_RK46_4, singlets_RK57_4, singlets_RK58_4, singlets_RK61_4)
length(cells.to.subset) #223285
object <- subset(object, cells = cells.to.subset)
object #25399 features across 223285 samples within 1 assay 



#================= Add metadata =================
# Add exp.cond
object@meta.data$exp.cond <- plyr::mapvalues(
  x = object@meta.data$orig.ident,
  from = c("RK1_4", "RK3_4", "RK5_4", "RK7_4", "RK9_4", "RK16_4", "RK43_4", "RK45_4", "RK46_4", "RK57_4", "RK58_4", "RK61_4"),
  to = c("Lean", "Lean", "Lean", "Obese", "Obese", "Obese", "Obese+sGCact", "Obese+sGCact", "Obese+sGCact", "Obese+sGCstim", "Obese+sGCstim", "Obese+sGCstim")
)

# Add uPCR
object@meta.data$uPCR <- plyr::mapvalues(
  x = object@meta.data$orig.ident,
  from = c("RK1_4", "RK3_4", "RK5_4", "RK7_4", "RK9_4", "RK16_4", "RK43_4", "RK45_4", "RK46_4", "RK57_4", "RK58_4", "RK61_4"),
  to = c("83", "85", "71", "3218", "2723", "2316", "445", "414", "395", "1043", "654", "1246")
)



#================= Save outfile =================
saveRDS(object, '/~/seurat/~/allcells_MASTER.rds')



#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)



