library(dplyr)
library(Seurat)
library(data.table)
library(cowplot)
library(ggplot2)
library(biomaRt)
library(Hmisc)
set.seed(123)



#============================================
#========= LOAD sGC rat REFERENCE ===========
#============================================
ZSF1rat <- readRDS('/~/seurat/~/object_processed_allcells.rds')



#========================================
#========= LOAD QUERY dataset ===========
#========================================
object <- readRDS("/~/GSE131882_Wilson_2019_humanearlyDKD/object.rds")

#add clinical/histopathology data
object@meta.data$eGFR <- plyr::mapvalues(
  x = object@meta.data$orig.ident,
  from = c("Control1", "Control2", "Control3",
           "DKD1", "DKD2", "DKD3"),
  to = c("58", "61", "69",
         "76", "56", "86")
)
object@meta.data$eGFR <- as.numeric(object@meta.data$eGFR)

object@meta.data$interst_fibrosis <- plyr::mapvalues(
  x = object@meta.data$orig.ident,
  from = c("Control1", "Control2", "Control3",
           "DKD1", "DKD2", "DKD3"),
  to = c("1-10%", "1-10%", "1-10%",
         "11-25%", "11-25%", "1-10%")
)

object@meta.data$glomerulosclerosis <- plyr::mapvalues(
  x = object@meta.data$orig.ident,
  from = c("Control1", "Control2", "Control3",
           "DKD1", "DKD2", "DKD3"),
  to = c("<10%", "<10%", "<10%",
         "10-25%", "26-50%", "<10%")
)

object@meta.data$exp.cond <- plyr::mapvalues(
  x = object@meta.data$orig.ident,
  from = c("Control1", "Control2", "Control3",
           "DKD1", "DKD2", "DKD3"),
  to = c("Control", "Control", "Control",
         "DKD", "DKD", "DKD")
)

#sort identities
object@meta.data$Idents <- factor(object@meta.data$Idents, 
                                  levels = c("Endo", "Podo","PEC","Mesangial", 
                                             "PCT", "LOH", "DCT", "DCT/CD", "PC", "IC_A", "IC_B", 
                                             "Leukocyte"))
Idents(object) <- object@meta.data$Idents
table(Idents(object))

#preprocessing
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
var.features <- object@assays$RNA@var.features
object <- ScaleData(object, features = var.features, vars.to.regress=c("nCount_RNA", "percent.mt"))
object <- RunPCA(object, features = VariableFeatures(object = object), seed.use=42)
object <- object %>% 
  RunUMAP(reduction = "pca", dims = 1:30, umap.method="uwot", n.neighbors = 30, n.components=2, metric="cosine", learning.rate=1, min.dist = 0.3,
          spread=1, set.op.mix.ratio=1, local.connectivity=1, repulsion.strength=1, negative.sample.rate=5, uwot.sgd=F, seed.use=42, angular.rp.forest=F) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30, k.param=20, compute.SNN=T, prune.SNN=0.0667, nn.method="rann", annoy.metric="euclidian", nn.eps=0,) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                              2.2, 2.4, 2.6, 2.8, 3.0,
                              4.0, 5.0, 6.0, 8.0, 10.0), save.SNN=TRUE) %>% 
  identity()



#===============================================================================
#===================== ORTHOLOGOUS MAPPING HUMAN TO RAT ========================
#===============================================================================
humangenes <- as.data.frame(rownames(object))
colnames(humangenes) <- "gene"
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                 filters = "hgnc_symbol", 
                 values = humangenes$gene , 
                 mart = human, 
                 attributesL = c("rgd_symbol"), 
                 martL = rat, 
                 uniqueRows=T)
genesV2
ratgenes <- genesV2$RGD.symbol
length(ratgenes) #13739

#subset object to those HUMAN genes for which a ratgene is available
GOI <- rownames(object)[rownames(object) %in% genesV2$HGNC.symbol]
length(GOI) #15450
rat_object <- subset(object, features = GOI)
rat_object #15450 features across 16259 samples within 1 assay 
target <- rownames(rat_object)
genesV2 <- genesV2[match(target, genesV2$HGNC.symbo), ]
counts <- GetAssayData(rat_object,assay = "RNA",slot = "counts")
rownames(counts) <- genesV2$Rat.gene.name
dim(counts) #15450 16259
`%notin%` <- Negate(`%in%`)
counts <- counts[rownames(counts) %notin% "",]
dim(counts) #12472 16259
WilsonDKD <- CreateSeuratObject(counts = counts,meta.data = rat_object@meta.data)
WilsonDKD #12472 features across 16259 samples within 1 assay
WilsonDKD@graphs <- rat_object@graphs
WilsonDKD@reductions <- rat_object@reductions
WilsonDKD@commands <- rat_object@commands
WilsonDKD <- NormalizeData(WilsonDKD, normalization.method = "LogNormalize", scale.factor = 10000)
var.features <- WilsonDKD@assays$RNA@var.features
WilsonDKD <- ScaleData(WilsonDKD, features = var.features, vars.to.regress=c("nCount_RNA", "percent.mt"))
#SUMMARY: WilsonDKD is a Seurat object with preprocessing, UMAP, clustering etc. as per the authors' instructions. 
#The cell-gene matrix was reduced to orthologous rat genes for downstream mapping purposes in Seurat.



#=============================================================
#========= MAPPING HUMAN DKD to ZSF1 DKD REFERENCE ===========
#=============================================================
#Find anchors
WilsonDKD <- WilsonDKD
anchors <- FindTransferAnchors(reference = ZSF1rat, 
                               query = WilsonDKD,
                               dims = 1:30,
                               reference.reduction = "pca") #Found 11293 anchors
saveRDS(anchors, 'anchors.rds')
predictions <- TransferData(anchorset = anchors, 
                            refdata = ZSF1rat$clusters2,
                            dims = 1:30)
saveRDS(predictions, 'predictions.rds')
WilsonDKD <- AddMetaData(WilsonDKD, metadata = predictions)

#Unimodal UMAP projection
ZSF1rat <- RunUMAP(ZSF1rat, dims = 1:30, reduction = "harmony", return.model = TRUE, n.neighbors=50, min.dist=0.9)
WilsonDKD <- MapQuery(anchorset = anchors, 
                      reference = ZSF1rat, 
                      query = WilsonDKD,
                      refdata = list(celltype = "clusters2"), 
                      reference.reduction = "harmony", 
                      reduction.model = "umap")
WilsonDKD@meta.data$predicted.celltype <- factor(WilsonDKD@meta.data$predicted.celltype, 
                                          levels = c("GEC(Meis2+)", "GEC(Reln+)", "AA/LA/EA", "Podo", 
                                                     "Fib/Peri", "SMC/Peri", 
                                                     "PCT", "PCT(Spp1+)", "PST(S2)", "PST", "PST(Spp1+)", "PTinj", "PTnov1",
                                                     "DTL", "DTL(Shroom3+)", "ATL", "TAL", 
                                                     "DCT", "DCT/CNT", "CNT/PC", 
                                                     "IC-A", "IC-A/CD-Trans", "IC-B", 
                                                     "Immune"))
Idents(WilsonDKD) <- WilsonDKD@meta.data$predicted.celltype
table(Idents(WilsonDKD))



#================================================
#========= CORRELATE CLUSTER AVERAGES ===========
#================================================
#load Wilson cluster averages
Idents(WilsonDKD) <- WilsonDKD$clusters1
cluster.averages <- AverageExpression(WilsonDKD, return.seurat = FALSE, 
                                      slot = "scale.data", 
                                      verbose = TRUE)
saveRDS(cluster.averages[["RNA"]], "/~/clusters1_cluster.averages_scaled__RATORTHOLOGUES.rds")
avg_Wilson <- as.data.frame(readRDS("/~/clusters1_cluster.averages_scaled__RATORTHOLOGUES.rds"))
avg_Wilson$X <- rownames(avg_Wilson)

#load sGC cluster averages
avg_sGC <- as.data.frame(readRDS("/~/ZSF1ratdownsample100k_MAST_clusters4_cluster.averages.rds"))
avg_sGC$X <- rownames(avg_sGC)
df_combined <- merge(avg_Wilson, avg_sGC, by = "X")
rownames(df_combined) <- df_combined$X
df_combined <- df_combined[ , !(names(df_combined) %in% c("X"))]
cor_PCC <- rcorr(as.matrix(df_combined))



#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)


