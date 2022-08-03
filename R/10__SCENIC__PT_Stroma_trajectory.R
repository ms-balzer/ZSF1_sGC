library(dplyr)
library(data.table)
library(cowplot)
library(ggplot2)
library(Seurat)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(doMC)
library(R2HTML)
library(rbokeh)
options(width=200)
set.seed(123)
dir.create("int")



### ========================== STEP A: Extract metadata from trajectory
PTstroma <- readRDS("/~/seurat/~/PTstroma.rds")
PTstroma@reductions$umap@assay.used <- "RNA"
PTstroma_traj <- subset(PTstroma, idents = c("PST", "PTinj", "ProfibPT", "DediffPT_1", "Int", "Mesench"), invert=FALSE, downsample=1027)
PTstroma_traj #25399 features across 4821 samples within 1 assay

#extract counts
counts <- as.sparse(PTstroma_traj@assays$RNA@counts)
saveRDS(counts, "int/exprMat.rds")

#extract cell annotation info
cellInfo <- data.frame(clusters2=PTstroma_traj$clusters2)
cellInfo$clusters2 <- as.factor(cellInfo$clusters2)
saveRDS(cellInfo, file="int/cellInfo.rds")

#extract embeddings
dr_coords <- PTstroma_traj@reductions$diffmap@cell.embeddings
saveRDS(dr_coords, file="int/dr_coords.rds")

#extract exp.cond
groupInfo <- data.frame(PTstroma_traj$exp.cond)
groupInfo$exp.cond <- as.factor(groupInfo$exp.cond)
saveRDS(groupInfo, file="int/groupInfo_exp.cond.rds")



### ========================== STEP B: Prep
exprMat <- as.matrix(counts)
org <- "mgi"
dbDir <- "/~/SCENIC/~" # RcisTarget databases location
myDatasetTitle <- "PT_Stroma_trajectory_GRN"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=12) 
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.rds"
scenicOptions@inputDatasetInfo$dr_coords <- "int/dr_coords.rds"
saveRDS(scenicOptions, file="int/scenicOptions.rds") 



### ========================== STEP C: Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
saveRDS(exprMat_filtered, 'int/exprMat_filtered.rds')
runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered, scenicOptions)



### ========================== STEP D: Build and score the GRN
scenicOptions <- readRDS("int/scenicOptions.rds")
genesKept <- loadInt(scenicOptions, "genesKept")
exprMat_filtered <- readRDS('int/exprMat_filtered.rds')
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 4
scenicOptions@settings$seed <- 123
exprMat_filtered <- as.matrix(exprMat_filtered)
exprMat_log <- log2(exprMat_filtered+1)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)



### ========================== STEP E: Get specific regulons
cellInfo <- readRDS("int/cellInfo.rds")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$clusters2),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), 1])
colnames(rss) <- c("Int", "PST", "ProfibPT", "Mesench", "PTinj_2", "DediffPT_1", "PTinj_1")
rssPlot <- plotRSS(rss, cluster_columns = FALSE)



#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)


