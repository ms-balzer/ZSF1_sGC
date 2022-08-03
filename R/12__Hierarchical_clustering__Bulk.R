library(ggplot2)
library(stringr)
library(dplyr)
library(broom)
library(data.table)
library(cluster)
library(factoextra)
library(tidyverse)
library(dendextend)
library(stringr)
library(ggpubr)
set.seed(123)



#======== STEP 1: load human data ========
meta <- read.csv("/data/Biobank-Tubule/tubule_metadata.csv", na.strings='.')
dim(meta) #991  43
dat <- read.csv("/~/data/Biobank-Tubule/HK.Biobank.Tubule.TPM.csv")
dim(dat) #44328   991



#======== STEP 2: get composite sGC co-expression WGCNA score genes and lift over from rat to human ========
geneInfo = read.csv('~/WGCNA/geneInfoSigned_UCI_all_by_clusters2.csv', sep = ",", header = TRUE)
geneInfo <- subset(geneInfo, Initially.Assigned.Module.Color%in%c("red", "green", "black","blue","yellow"))
genelist <- geneInfo$X
length(unique(genelist)) #2198
ratgenes <- as.data.frame(unique(genelist))
colnames(ratgenes) <- "gene"
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
genesV2 = getLDS(attributes = c("rgd_symbol"), 
                 filters = "rgd_symbol", 
                 values = ratgenes$gene , 
                 mart = rat, 
                 attributesL = c("hgnc_symbol"), 
                 martL = human, 
                 uniqueRows=T)
genesV2
humangenes <- unique(genesV2$HGNC.symbol)
length(humangenes) #1901



#======== STEP 3: do cluster analysis ========
#subset TPM matrix on humangenes
dat <- dat[which(rownames(dat)%in% humangenes),]
df <- as.data.frame(dat)
dim(df) #1755  991
df2 <- t(df)
df_sc <- as.data.frame(scale(df2))

#create distance matrix and dendrogram
dist_mat <- dist(df_sc, method ="euclidean")
hclust_avg <- hclust(dist_mat, method="ward")
pdf('dendrogram.pdf')
plot(hclust_avg)
dev.off()

#determine optimal number of clusters and cut dendrogram
pdf('optimal_n_of_clusters.pdf')
fviz_nbclust(df2, kmeans, method = "silhouette")
dev.off()
k=2 #choose k based on plot above
cut_avg <- cutree(hclust_avg, k=k)
pdf(paste0('dendrogram_with_clustering_k',k,'_color.pdf'), width=6, height=4)
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = k, col = c("#990000", "gray"))
plot(avg_col_dend)
dev.off()



#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)


