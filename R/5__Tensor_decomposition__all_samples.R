library(scITD)
library(Seurat)
library(dplyr)
library(data.table)
library(RColorBrewer)
options(scipen=1000000)
set.seed(123)



#================ LOAD DATA ================
object <- readRDS('/~/seurat/~/object_processed_allcells.rds')

# counts
counts <- readRDS('counts.rds')

# meta data
meta <- object@meta.data
meta$ctypes = meta$clusters3
meta$donors = meta$orig.ident

numbers_of_bins = 3
meta<-meta%>%mutate(uPCR_bins = cut(uPCR, 
                                    breaks = unique(quantile(uPCR,probs=seq.int(0,1, by=1/numbers_of_bins))), 
                                    include.lowest=TRUE))
meta$uPCR_bins <- plyr::mapvalues(x = meta$uPCR_bins,
                                  from = names(table(meta$uPCR_bins)),
                                  to = c("71-414", "415-1250", "1251-3220"))
keep <- c("donors","ctypes",
          "interst_fibrosis", "tub_degen","hyaline_cast","mononuc_infiltr","glomerulopathy","uPCR_bins","genotype","treatment")
meta = meta[keep]
saveRDS(meta, 'meta.rds')



#================ SETUP ================
# set up project parameters
param_list <- initialize_params(ctypes_use = c("Endo",
                                               "Podo", 
                                               "Stroma",
                                               "Prox tub",
                                               "LOH",
                                               "DCT/CNT/PC",
                                               "IC",
                                               "Immune"),
                                ncores = 8, rand_seed = 10)

# create project container
container <- make_new_container(count_data=counts, 
                                     meta_data=meta,
                                     params=param_list,
                                     label_donor_sex = F)

# form the tensor from the data
container <- form_tensor(container, 
                              donor_min_cells=5,
                              norm_method='trim', 
                              scale_factor=10000,
                              vargenes_method='norm_var_pvals',
                              vargenes_thresh=.1,
                              scale_var = TRUE, 
                              var_scale_power = 2)

# number of genes included in the tensor
#check the number of overdispersed genes identified to make sure the tensor has a decent number of genes before running the decomposition
#To increase the number of genes included, simply rerun form_tensor() with a higher vargenes_thresh value.
print(length(container[["all_vargenes"]])) #1430

# get assistance with rank determination
container <- determine_ranks_tucker(container, 
                                         max_ranks_test=c(10,15),
                                         shuffle_level='cells', 
                                         num_iter=10, 
                                         norm_method='trim',
                                         scale_factor=10000,
                                         scale_var=TRUE,
                                         var_scale_power=2)
pdf('rank_determination.pdf')
container$plots$rank_determination_plot
dev.off()

# run stability analysis
factors <- 4
container <- run_stability_analysis(container,
                                         ranks=c(factors,10),
                                         n_iterations=50,
                                         subset_type='subset', 
                                         sub_prop=.95)
pdf(paste0('stability_analysis_',factors,'_factors.pdf'))
container$plots$stability_plot_dsc
dev.off()



#================ RUN TUCKER TENSOR DECOMPOSITION ================
# run the tensor decomposition
container <- run_tucker_ica(container, 
                            ranks=c(factors,10),
                            tucker_type = 'regular', 
                            rotation_type = 'hybrid')



#================ PLOT DONOR SCORES MATRIX (WHICH FACTOR IS PRESENT IN EACH DONOR?) ================
# get donor scores-metadata associations
container <- get_meta_associations(container,
                                   vars_test=c("interst_fibrosis", "tub_degen","hyaline_cast","mononuc_infiltr","glomerulopathy","uPCR_bins","genotype","treatment"),
                                   stat_use='pval')

# plot donor scores
container <- plot_donor_matrix(container,
                               meta_vars=c("interst_fibrosis", "tub_degen","hyaline_cast","mononuc_infiltr","glomerulopathy","uPCR_bins","genotype","treatment"),
                               cluster_by_meta = 'interst_fibrosis',
                               show_donor_ids = TRUE,
                               add_meta_associations='pval')

container$plots$donor_matrix@ht_list$interst_fibrosis@matrix_color_mapping@colors <- brewer.pal(n = 5, name = "YlOrRd")
names(container$plots$donor_matrix@ht_list$interst_fibrosis@matrix_color_mapping@colors) <- c("0", "1", "2", "3", "4")

dummypalette <- brewer.pal(n = 4, name = "Oranges")
container$plots$donor_matrix@ht_list$tub_degen@matrix_color_mapping@colors <- c(dummypalette[2], dummypalette[1], dummypalette[3:4])
names(container$plots$donor_matrix@ht_list$tub_degen@matrix_color_mapping@colors) <- c("1", "0", "2", "3")

container$plots$donor_matrix@ht_list$hyaline_cast@matrix_color_mapping@colors <- brewer.pal(n = 4, name = "PuRd")
names(container$plots$donor_matrix@ht_list$hyaline_cast@matrix_color_mapping@colors) <- c("0", "1", "2", "3")

container$plots$donor_matrix@ht_list$mononuc_infiltr@matrix_color_mapping@colors <- brewer.pal(n = 4, name = "Greens")
names(container$plots$donor_matrix@ht_list$mononuc_infiltr@matrix_color_mapping@colors) <- c("0", "1", "2", "3")

container$plots$donor_matrix@ht_list$glomerulopathy@matrix_color_mapping@colors <- brewer.pal(n = 5, name = "YlGnBu")
names(container$plots$donor_matrix@ht_list$glomerulopathy@matrix_color_mapping@colors) <- c("0", "1", "2", "3", "4")

container$plots$donor_matrix@ht_list$uPCR_bins@matrix_color_mapping@colors <- brewer.pal(n = 3, name = "YlGnBu")
names(container$plots$donor_matrix@ht_list$uPCR_bins@matrix_color_mapping@colors) <- c("71-414", "415-1250", "1251-3220")

container$plots$donor_matrix@ht_list$genotype@matrix_color_mapping@colors <- c("#E0E0E0", "#000000")
names(container$plots$donor_matrix@ht_list$genotype@matrix_color_mapping@colors) <- c("Lean", "Obese")

container$plots$donor_matrix@ht_list$treatment@matrix_color_mapping@colors <- c("Orange", "Purple")
names(container$plots$donor_matrix@ht_list$treatment@matrix_color_mapping@colors) <- c("No sGCm", "sGCm treatment")

# show the donor scores heatmap
pdf(paste0('donor_scores_heatmap_',factors,'_factors.pdf'), width=8)
container$plots$donor_matrix
dev.off()

# show how much each factor explains in variation
container$exp_var
# 48.710242 18.959396  8.263184  4.141697



#================ DETERMINE WHICH GENES FROM EACH CELL TYPE ARE SIGNIFICANTLY ASSOCIATED W/ EACH FACTOR ================
# get significant genes
container <- get_lm_pvals(container)

#extract the p values for the gene-factor associations
pvalues <- container$gene_score_associations 
pvalues <- as.data.frame(pvalues)
write.table(x = pvalues, file = paste0("gene_score_asccociations_",factors,"_factors.csv"), 
            append = FALSE, sep = ",", dec = ".", row.names = TRUE, col.names = NA)

# generate the loadings plots
container <- get_all_lds_factor_plots(container, 
                                      use_sig_only=F,
                                      nonsig_to_zero=F,
                                      annot = 'sig_genes',
                                      sig_thresh=.05,
                                      display_genes=F,
                                      gene_callouts = TRUE,
                                      callout_n_gene_per_ctype=5,
                                      show_var_explained = TRUE,
                                      reset_other_factor_plots = TRUE)



#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)


