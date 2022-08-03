library(slingshot)
library(Seurat)
library(ggbeeswarm)
library(ggthemes)
library(SingleCellExperiment)
library(RColorBrewer)
library(destiny)
library(gam)
set.seed(123)



#========================================================================================
#================== LOAD PT & STROMA, SUBSET TO TRAJECTORY, CALC DIM RED ================
#========================================================================================
PTstroma <- readRDS("/~/seurat/~/PTstroma.rds")
PTstroma@reductions$umap@assay.used <- "RNA"
PTstroma_traj <- subset(PTstroma, idents = c("PST", "PTinj", "ProfibPT", "DediffPT_1", "Int", "Mesench"), invert=FALSE, downsample=1027)
PTstroma_traj #25399 features across 4821 samples within 1 assay

pcs = PTstroma_traj@reductions$pca
emb.pca = pcs@cell.embeddings
diffmap = destiny::DiffusionMap(emb.pca, k=50)
emb.diffmap = destiny::eigenvectors(diffmap)[,1:2]
rownames(emb.diffmap) <- rownames(emb.pca)

PTstroma_traj@reductions$diffmap <- PTstroma_traj@reductions$umap
PTstroma_traj@reductions$diffmap@cell.embeddings <- emb.diffmap
PTstroma_traj@reductions$diffmap@key <- "DC"




#==========================================================================
#================== CONVERT SEURAT TO SINGLECELLEXPERIMENT ================
#==========================================================================
sce <- as.SingleCellExperiment(PTstroma_traj)
sce@colData$orig.ident <- factor(sce@colData$orig.ident, levels = c("RK1_4", "RK3_4", "RK5_4", 
                                                                    "RK7_4", "RK9_4", "RK16_4",
                                                                    "RK43_4", "RK45_4", "RK46_4", 
                                                                    "RK57_4", "RK58_4", "RK61_4"))
sce@colData$exp.cond <- factor(sce@colData$exp.cond, levels = c("Lean", "Obese","Obese+sGCact","Obese+sGCstim"))
sce@colData$clusters1 <- factor(sce@colData$clusters1, levels = c("PST", "PTinj", "ProfibPT", "DediffPT_1", "Int", "Mesench"))



#===============================================================
#================== SLINGSHOT UPSTREAM ANALYSIS ================
#===============================================================

# ==== Gene filtering ====
# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]



# ==== Normalization ====
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)



# ==== Identify global lineage structure ====
lin0c <- getLineages(sce@int_colData@listData$reducedDims@listData$DIFFMAP, sce$clusters1, start.clus= c('ProfibPT'))
lin0c



# ==== Construct smooth curves and ordering cells ====
crv0c <- getCurves(lin0c)
crv0c #inspect, need to edit the curves and get rid of the weird ones
crv0cedit <- crv0c
crv0cedit@lineages$Lineage2 <- NULL
crv0cedit@curves$curve2 <- NULL
crv0cedit



# ==== Extract pseudotime from crv ====
pt0cedit <- slingPseudotime(crv0cedit, na=TRUE)
summary(pt0cedit)

#PTinj splits into 2 separate clusters on lineages 1 and 2, respectively. Rename "PTinj" cells in lineages 1 and 2 as "PTinj_1" and "PTinj_2", respectively.
testdf <- as.data.frame(pt0cedit) #convert to df
testdf$mean <- rowMeans(testdf[,c('curve1', 'curve3')], na.rm=TRUE) #some cells are in both lineages, so we need to take the mean of the pseudotime value
#store in sce object
sce@colData$pt0cedit.1 <- pt0cedit[,1]
sce@colData$pt0cedit.2 <- pt0cedit[,2]
sce@colData$pt0cedit_common <- testdf$mean
names(sce@colData$pt0cedit_common) <- names(sce@colData$pt0cedit.2)
#get lin2 PTinj barcodes (former curve 1)
rownames(testdf[!is.na(testdf$curve1),]) #barcodes in curve1, n=2913
rownames(testdf[sce$clusters1=="PTinj",]) #barcodes for PTinj n=1019
crv1_PTinj <- rownames(testdf[sce$clusters1=="PTinj" & !is.na(testdf$curve1),]) #barcodes for both conditions, n=541
table(sce$clusters1)
#PST      PTinj   ProfibPT DediffPT_1        Int         MC 
#917       1019       1027        534        306       1018
#export factor to df
dummy <- as.character(sce$clusters1) #create character vector of clusters1
names(dummy) <- colnames(sce) #name with barcodes
dummy[sce$clusters1=="PTinj" & !is.na(testdf$curve1)] <- "PTinj_2" #rename those in question
sce$clusters2 <- as.factor(dummy) #copy back as factor
sce$clusters2 <- plyr::mapvalues(
  x = sce$clusters2,
  from = c("DediffPT_1", "Int", "Mesench", "ProfibPT", "PST", "PTinj", "PTinj_2"),
  to = c("DediffPT_1", "Int", "Mesench", "ProfibPT", "PST", "PTinj_1", "PTinj_2"))# recode PTinj as PTinj_1
sce$clusters2 <- factor(sce$clusters2, levels = c("PST", "PTinj_1", "ProfibPT", "PTinj_2", "DediffPT_1", "Int", "Mesench"))



#=================================================================
#================== SLINGSHOT DOWNSTREAM ANALYSIS ================
#=================================================================

# ==== Identify temporally expressed genes (using GAM) ====

t <- sce$pt0cedit.1[!is.na(sce$pt0cedit.1)]
sce_lin2 <- sce[,names(which(!is.na(sce$pt0cedit.1))==TRUE)]
Y <- log1p(assays(sce_lin2)$norm)
var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y <- Y[var100,]
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(sce_lin2)$norm[topgenes, order(t, na.last = NA)]
heatclus <- sce_lin2$clusters2[order(t, na.last = NA)]

# ==== Create pt bins
df<-data.frame(pseudotime = sce_lin2$pt0cedit.1)
numbers_of_bins = 6
df<-df%>%mutate(MyQuantileBins = cut(pseudotime, 
                                     breaks = unique(quantile(pseudotime,probs=seq.int(0,1, by=1/numbers_of_bins))), 
                                     include.lowest=TRUE))
df$MyQuantileBins <- plyr::mapvalues(
  x = df$MyQuantileBins,
  from = c(names(table(df$MyQuantileBins))[1], names(table(df$MyQuantileBins))[2], names(table(df$MyQuantileBins))[3], 
           names(table(df$MyQuantileBins))[4], names(table(df$MyQuantileBins))[5], names(table(df$MyQuantileBins))[6]),
  to = c("pt1", "pt2", "pt3", "pt4", "pt5", "pt6"))
sce_lin2$pt_lin2_bin <- df$MyQuantileBins



#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)


