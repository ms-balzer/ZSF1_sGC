#This pipeline is run individually on every sample, here we use sample RK1_4_comb.
library(dplyr)
library(Seurat)
library(SoupX)
library(data.table)
library(DoubletFinder)
set.seed(123)



#===============================
#============ SoupX ============
#===============================
RK1_4_comb.data <- Read10X(data.dir = "/~/cellranger_count/~/outs/filtered_feature_bc_matrix")
RK1_4_comb <- CreateSeuratObject(counts = RK1_4_comb.data, project = "RK1_4_comb")
sc_RK1_4_comb = load10X("/~/cellranger_count/~/outs/")
sc_RK1_4_comb = autoEstCont(sc_RK1_4_comb)
out_RK1_4_comb = adjustCounts(sc_RK1_4_comb, roundToInt=T)
RK1_4_comb@assays$RNA@counts@x <- out_RK1_4_comb@x
RK1_4_comb@assays$RNA@data@x <- out_RK1_4_comb@x



#=========================================
#============ Seurat pipeline ============
#=========================================
RK1_4_comb[["percent.mt"]] <- PercentageFeatureSet(RK1_4_comb, pattern = "^Mt-")
RK1_4_comb <- subset(RK1_4_comb, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15)
RK1_4_comb <- NormalizeData(RK1_4_comb, normalization.method = "LogNormalize", scale.factor = 10000)
RK1_4_comb <- FindVariableFeatures(RK1_4_comb, selection.method = "vst", nfeatures = 3000)
RK1_4_comb <- ScaleData(RK1_4_comb, features = rownames(RK1_4_comb))
RK1_4_comb <- RunPCA(RK1_4_comb, features = VariableFeatures(object = RK1_4_comb))
RK1_4_comb <- RK1_4_comb %>% 
  RunUMAP(reduction = "pca", dims = 1:15) %>% 
  FindNeighbors(reduction = "pca", dims = 1:15) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0), save.SNN=TRUE) %>% 
  identity()



#=========================================
#============ Remove Doublets ============
#=========================================
# pK Identification
sweep.res.list_kidney <- paramSweep_v3(RK1_4_comb, PCs = 1:15, sct = FALSE, num.cores=8)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
pK <- bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]
pK <- as.numeric(as.vector(pK))

# Homotypic Doublet Proportion Estimate
annotations <- RK1_4_comb@meta.data$RNA_snn_res.1.1
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(homotypic.prop*length(colnames(RK1_4_comb)))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies
RK1_4_comb <- doubletFinder_v3(RK1_4_comb, PCs = 1:15, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
slotname <- paste0("DF.classifications_0.25_",pK,"_",nExp_poi.adj)
table(RK1_4_comb$"DF.classifications_0.25_0.001_982")
saveRDS(RK1_4_comb, '~/seurat/RK1_4_comb_soupX/object_doubletfinder.rds')



#========= Subtract doublets and plot again ===========
Idents(RK1_4_comb) <- RK1_4_comb$"DF.classifications_0.25_0.001_982"
singlets <- WhichCells(RK1_4_comb, idents = "Singlet")
doublets <- WhichCells(RK1_4_comb, idents = "Doublet")
saveRDS(singlets, '/~/seurat/RK1_4_comb_soupX/singlets.rds')
saveRDS(doublets, '/~/seurat/RK1_4_comb_soupX/doublets.rds')



#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)


