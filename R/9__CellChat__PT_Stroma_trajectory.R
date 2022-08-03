library(CellChat)
library(patchwork)
library(data.table)
library(Seurat)
library(ggplot2)
library(NMF)
library(ggalluvial)
library(viridis)
library(reticulate)
options(stringsAsFactors = FALSE)
set.seed(123)



#==========================================================================================
#==========================================================================================
#========= Part I: Data input & processing and initialization of CellChat object ==========
#==========================================================================================
#==========================================================================================

#=============== Load data
PTstroma <- readRDS("/~/seurat/~/PTstroma.rds")
PTstroma@reductions$umap@assay.used <- "RNA"
PTstroma_traj <- subset(PTstroma, idents = c("PST", "PTinj", "ProfibPT", "DediffPT_1", "Int", "Mesench"), invert=FALSE, downsample=1027)
PTstroma_traj #25399 features across 4821 samples within 1 assay
col.palette <- c("#336600", "#FF0000", "#660000", "#990099", "#FF66FF", "#404040", "#A0A0A0")



#=============== Create a CellChat object
cellchat <- createCellChat(object = PTstroma_traj, group.by = "clusters2", assay = "RNA")
cellchat #An object of class CellChat created from a single dataset 
#25399 genes.
#4821 cells.



#=============== Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use



#=============== Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 8) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI network
cellchat <- projectData(cellchat, PPI.mouse)





#==========================================================================
#==========================================================================
#========= Part II: Inference of cell-cell communication network ==========
#==========================================================================
#==========================================================================

#=============== Compute the communication probability and infer cellular communication network
#The function computeAveExpr can help to check the average expression of signaling genes of interest, e.g, 
cellchat <- computeCommunProb(cellchat,
                              type='triMean',
                              trim = NULL)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)



#=============== Extract the inferred cellular communication network as a data frame
#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
                                        #Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
fwrite(x = df.net, row.names = TRUE, file = 'df.net.csv')



#=============== Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)



#=============== Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat) #can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use





#===============================================================================
#===============================================================================
#========= Part III: Visualization of cell-cell communication network ==========
#===============================================================================
#===============================================================================

#=============== Automatically save the plots of the all inferred network for quick exploration
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways

# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,7)

#for all significant pathways
pathways.show.all <- cellchat@netP$pathways





#=================================================================================
#=================================================================================
#========= Part IV: Systems analysis of cell-cell communication network ==========
#=================================================================================
#=================================================================================

#=======================================================================================================================================
#=============== Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
#=======================================================================================================================================

#=============== Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pathways.show = pathways.show.all
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10, color.use=col.palette)



#=============== Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat, color.use=col.palette, label.size = 4)



#=============== Identify signals contributing most to outgoing or incoming signaling of certain cell groups
#We can also answer the question on which signals contributing most to outgoing or incoming signaling of certain cell groups.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.use=col.palette, color.heatmap = "Purples", font.size = 12, font.size.title=12, width=5)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.use=col.palette, color.heatmap = "Purples", font.size = 12, font.size.title=12, width=5)



#====================================================================================================================================
#=============== Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
#====================================================================================================================================

#=============== Identify and visualize outgoing communication pattern of secreting cells
#Here we run selectK to infer the number of patterns.
selectK(cellchat, pattern = "outgoing")
nPatterns = 2 #pick based on k value
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, color.heatmap = "PuOr", color.use=col.palette, font.size = 12, width=2, height=5)
netAnalysis_river(cellchat, pattern = "outgoing", color.use=col.palette, 
                  color.use.pattern = c("Pattern 1"="red", "Pattern 2"="yellow", "Pattern 3"="blue"), 
                  color.use.signaling=viridis(14),font.size = 2.5)
netAnalysis_dot(cellchat, pattern = "outgoing", color.use=col.palette, dot.size = c(2, 6),font.size = 14,font.size.title = 14)



#=============== Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, color.heatmap = "PuOr", color.use=col.palette[c(1:5,7)],font.size = 12, width=2, height=5)
netAnalysis_river(cellchat, pattern = "incoming", color.use=col.palette, 
                  color.use.pattern = c("Pattern 1"="red", "Pattern 2"="yellow", "Pattern 3"="blue"), 
                  color.use.signaling=viridis(14))
netAnalysis_dot(cellchat, pattern = "incoming", color.use=col.palette, dot.size = c(2, 6),font.size = 14,font.size.title = 14)



#===================================================================================
#=============== Manifold and classification learning analysis of signaling networks
#===================================================================================

#=============== Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional") # creates 'estimationNumCluster__functional_dataset_single.pdf'
netVisual_embedding(cellchat, type = "functional", title = "Functional similarity", label.size = 3.5, dot.size = c(2, 6), font.size = 10, font.size.title = 10, show.axes = T, show.legend = T)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

#=============== Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural") #creates 'estimationNumCluster__structural_dataset_single.pdf'
netVisual_embedding(cellchat, type = "structural", title = "Structural similarity", label.size = 3.5, dot.size = c(2, 6), font.size = 10, font.size.title = 10, show.axes = T, show.legend = T)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)





#=====================================================
#=====================================================
#========= Part V: Save the CellChat object ==========
#=====================================================
#=====================================================
saveRDS(cellchat, file = "cellchat_object.rds")





#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)


