library(dplyr)
library(Seurat)
library(data.table)
library(cowplot)
library(ggplot2)
library(harmony)
set.seed(123)



#========================================
#============ START ANALYSIS ============
#========================================

#================= Standard pre-processing workflow =================
object <- readRDS('/~/seurat/~/allcells_MASTER.rds')
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^Mt-")
object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(object)
object <- ScaleData(object, features = all.genes, vars.to.regress="nCount_RNA")
object <- RunPCA(object, features = VariableFeatures(object = object))



#============= Run Harmony =============
object <- object %>% RunHarmony("orig.ident", plot_convergence = TRUE)



#============= Downstream analysis =============
object <- object %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, n.neighbors = 50L, min.dist = 0.9) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4, 1.5), save.SNN=TRUE) %>% 
  identity()
object@meta.data$orig.ident <- factor(object@meta.data$orig.ident, 
                                           levels = c("RK1_4", "RK3_4", "RK5_4",
                                                      "RK7_4", "RK9_4", "RK16_4",
                                                      "RK43_4", "RK45_4", "RK46_4",
                                                      "RK57_4", "RK58_4", "RK61_4"))
Idents(object) <- object@meta.data$orig.ident
table(Idents(object))



#=================================================
#========= PROCEED WITH RESOLUTION 0.9 ===========
#=================================================
#clusters2
object@meta.data$clusters2 <- plyr::mapvalues(
  x = object@meta.data$RNA_snn_res.0.9,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 
           "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
           "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
           "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", 
           "40"),
  to = c("PST", "PCT", "PST(Spp1+)", "TAL", "GEC(Meis2+)", "PST(S2)", "CNT/PC", "DCT", "PTnov1", "PST", 
         "PST", "PCT(Spp1+)", "PST", "PTinj", "IC-A", "IC-B", "DCT/CNT", "Fib/Peri", "TAL", "PCT",
         "DTL(Shroom3+)", "Immune", "Podo", "PST", "PCT", "DTL", "TAL", "AA/LA/EA", "SMC/Peri", "PTnov2",
         "ATL", "Immune", "DCT", "GEC(Reln+)", "IC-A/CD-Trans", "DCT/CNT", "PCT", "PST", "Podo", "TAL", 
         "PTnov2"),
)
object@meta.data$clusters2 <- factor(object@meta.data$clusters2, 
                                          levels = c("GEC(Meis2+)", "GEC(Reln+)", "AA/LA/EA", "Podo", 
                                                     "Fib/Peri", "SMC/Peri", 
                                                     "PCT", "PCT(Spp1+)", "PST(S2)", "PST", "PST(Spp1+)", "PTinj", "PTnov1", "PTnov2", 
                                                     "DTL", "DTL(Shroom3+)", "ATL", "TAL", 
                                                     "DCT", "DCT/CNT", "CNT/PC", 
                                                     "IC-A", "IC-A/CD-Trans", "IC-B", 
                                                     "Immune"))



#================= Finding differentially expressed features (clusters2) =================
#find markers for every cluster compared to all remaining cells
object@meta.data$clusters2 <- factor(object@meta.data$clusters2, 
                                          levels = c("GEC(Meis2+)", "GEC(Reln+)", "AA/LA/EA", "Podo", 
                                                     "Fib/Peri", "SMC/Peri", 
                                                     "PCT", "PCT(Spp1+)", "PST(S2)", "PST", "PST(Spp1+)", "PTinj", "PTnov1", "PTnov2", 
                                                     "DTL", "DTL(Shroom3+)", "ATL", "TAL", 
                                                     "DCT", "DCT/CNT", "CNT/PC", 
                                                     "IC-A", "IC-A/CD-Trans", "IC-B", 
                                                     "Immune"))
Idents(object) <- object@meta.data$clusters2
table(Idents(object))

object.markers <- FindAllMarkers(object, test.use="MAST", only.pos = F, min.pct = 0.1, logfc.threshold = 0.2)
object.markers %>% group_by(cluster) -> tibble_FindAllMarkers_MAST



#===============================================================
#========= Add coarse-grained clustering (clusters3) ===========
#===============================================================
#add more coarse-grained clustering
object@meta.data$clusters3 <- plyr::mapvalues(
  x = object@meta.data$clusters2,
  from = c("GEC(Meis2+)", "GEC(Reln+)", "AA/LA/EA", "Podo", 
           "Fib/Peri", "SMC/Peri", 
           "PCT", "PCT(Spp1+)", "PST(S2)", "PST", "PST(Spp1+)", "PTinj", "PTnov1", "PTnov2", 
           "DTL", "DTL(Shroom3+)", "ATL", "TAL", 
           "DCT", "DCT/CNT", "CNT/PC", 
           "IC-A", "IC-A/CD-Trans", "IC-B", 
           "Immune"),
  to = c("Endo", "Endo", "Endo", 
         "Podo", 
         "Stroma", "Stroma", 
         "Prox tub", "Prox tub", "Prox tub", "Prox tub", "Prox tub", "Prox tub", "Prox tub", "Prox tub", 
         "LOH", "LOH", "LOH", "LOH", 
         "DCT/CNT/PC", "DCT/CNT/PC", "DCT/CNT/PC", 
         "IC", "IC", "IC", 
         "Immune")
)

#save outfile
saveRDS(object, '/~/seurat/~/object_processed_allcells.rds')



#============================================================================
#========= Export cell types for PT/Stroma subclustering analysis ===========
#============================================================================
Idents(object) <- object$clusters2
table(Idents(object))
object_PTstroma <- subset(object, idents = c("PCT", "PCT(Spp1+)", "PST(S2)", "PST", "PST(Spp1+)", "PTinj", "PTnov1", "PTnov2",
                                               "Fib/Peri", "SMC/Peri", "DTL(Shroom3+)"), invert = FALSE)
saveRDS(object_PTstroma, "/~/seurat/~/PTstroma.rds")



#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)


