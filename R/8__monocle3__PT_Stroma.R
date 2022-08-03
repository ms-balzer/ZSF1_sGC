library(dplyr)
library(Seurat)
library(data.table)
library(cowplot)
library(ggplot2)
library(monocle3)
library(htmlwidgets)
set.seed(123)



#=========================================================================
#========= LOADING DATA, CREATING CDS & ADDING METADATA ==================
#=========================================================================
### load cells
PTstroma <- readRDS("/~/seurat/~/PTstroma.rds")
PTstroma@reductions$umap@assay.used <- "RNA"
PTstroma_traj <- subset(PTstroma, idents = c("PST", "PTinj", "ProfibPT", "DediffPT_1", "Int", "Mesench"), invert=FALSE, downsample=1027)
PTstroma_traj #25399 features across 4821 samples within 1 assay



#====================================================================
#================= Normalize & pre-process the data =================
#====================================================================
cds <- preprocess_cds(cds, num_dim = 100, cores=8)
cds <- reduce_dimension(cds, 
                        umap.fast_sgd = FALSE, 
                        umap.min_dist = 0.8,
                        umap.n_neighbors = 300L,
                        cores=1)



#=====================================================
#================= Cluster the cells =================
#=====================================================
res=3e-4
cds = cluster_cells(cds, resolution=res)



#=========================================================================
#================= Constructing single-cell trajectories =================
#=========================================================================
cds <- learn_graph(cds, use_partition = FALSE)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, assigned_cell_type="PST"){
  cell_ids <- which(colData(cds)[, "assigned_cell_type"] == assigned_cell_type)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))



#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)


