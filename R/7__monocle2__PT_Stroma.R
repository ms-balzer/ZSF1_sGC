library(dplyr)
library(Seurat)
library(data.table)
library(cowplot)
library(ggplot2)
library(monocle)
library(htmlwidgets)
library(OneR)
set.seed(123)



#=========================================================================
#========= LOADING DATA, CREATING CDS & ADDING METADATA ==================
#=========================================================================
### load cells
PTstroma <- readRDS("/~/seurat/~/PTstroma.rds")
PTstroma@reductions$umap@assay.used <- "RNA"
PTstroma_traj <- subset(PTstroma, idents = c("PST", "PTinj", "ProfibPT", "DediffPT_1", "Int", "Mesench"), invert=FALSE, downsample=1027)
PTstroma_traj #25399 features across 4821 samples within 1 assay

### create CDS object
gene_annotation <- as.data.frame(PTstroma_traj@assays[["RNA"]]@counts@Dimnames[[1]], 
                                 row.names = PTstroma_traj@assays[["RNA"]]@counts@Dimnames[[1]])
gene_annotation <- AnnotatedDataFrame(gene_annotation)
colnames(gene_annotation) <- "gene_short_name"
cell_metadata <- as.data.frame(PTstroma_traj@assays[["RNA"]]@counts@Dimnames[[2]], 
                               row.names = PTstroma_traj@assays[["RNA"]]@counts@Dimnames[[2]])
cell_metadata <- AnnotatedDataFrame(cell_metadata)
colnames(cell_metadata) <- "barcode"
New_matrix <- PTstroma_traj@assays[["RNA"]]@counts
expression_matrix <- New_matrix
cds <- newCellDataSet(as(expression_matrix, "sparseMatrix"),
                                  phenoData = cell_metadata,
                                  featureData = gene_annotation,
                                  expressionFamily=negbinomial.size())

### add metadata
load(file="/home~/seurat/~/EXPORTtraj.RData")
cds@phenoData@data$percent.mt <- export__percent.mt
cds@phenoData@data$uPCR <- export__uPCR
cds@phenoData@data$orig.ident <- export__orig.ident
cds@phenoData@data$exp.cond <- export__exp.cond
cds@phenoData@data$clusters1 <- export__clusters1



#=================================================
#================= PREPROCESSING =================
#=================================================
#Estimate size factors and dispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#Filtering low-quality cells
cds <- detectGenes(cds, min_expr = 0.1)
cds #still 4821 cells
print(head(fData(cds)))



#==================================================================
#================= Classifying and Counting Cells =================
#==================================================================

#======= Clustering cells without marker genes 
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
cds <- reduceDimension(cds, max_components = 2, num_dim = 15,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds)



#=========================================================================
#================= Constructing Single Cell Trajectories =================
#=========================================================================
cds <- reduceDimension(cds, 
                       max_components = 2,
                       reduction_method = 'DDRTree', 
                       verbose = T)
cds <- orderCells(cds, reverse = F)



#====================================================================
#================= Differential Expression Analysis =================
#====================================================================
# Finding Genes that Distinguish Cell State
marker_genes <- row.names(fData(cds))
diff_test_res_State <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~State")
diff_test_res_State[,c("gene_short_name", "pval", "qval")]
sig_gene_names_State <- row.names(subset(diff_test_res_State, qval < 0.01))



# Clustering Genes by Pseudotemporal Expression Pattern
marker_genes <- row.names(fData(cds))
diff_test_res_Pseudotime <- differentialGeneTest(cds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names_Pseudotime <- row.names(subset(diff_test_res_Pseudotime, qval < 0.01))



#=========================================================================
#================= Analyzing Branches in Single-Cell Trajectories ========
#=========================================================================
BEAM_res <- BEAM(cds, branch_point = 3, cores = 8)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]



#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)


