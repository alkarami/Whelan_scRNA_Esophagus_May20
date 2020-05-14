## Load monocle and Seurat.
library(monocle)
library(Seurat)
library(cowplot)

## Load epithelial data.
load('ovy.epi')

## Use the integrated assay for pseudotime, and subset out non-basal cells.
ovy.b <- subset(ovy.epi, subset = c('Basal-1','Basal-2','Basal-3','Basal-4','Basal-5'))

## Using the monocon function defined in Monocle_PT_Epi, convert the basal cells for pseudotime with Monocle
bmono <- monocon(ovy.b)

## Visualize trajectory and the localization of proliferation marker Mki67. 
p1 <- plot_cell_trajectory(bmono, color_by = 'NamedTypes')
p2 <- plot_cell_trajectory(bmono, markers = 'Mki67', use_color_gradient = T)
plot_grid(p1,p2)

## Root trajectory at the localization of proliferation (State 10)
bmono <- orderCells(bmono, root_state = 10)

## 5 branchpoints should be identified. Perform BEAM 

############################### BEAMer performs BEAM and returns significantly branch-dependent genes (p<0.05) #########################

BEAMer <- function(object,branchpoint){
  BEAM_x <- BEAM(object, branch_point = branchpoint, cores = 4)
  BEAM_x <- BEAM_x[order(BEAM_x$qval),]
  BEAM_x <- BEAM_x[,c("gene_short_name", "pval", "qval")]
  BEAM_x_sig <- BEAM_x[which(BEAM_x$qval<0.05),]
  return(BEAM_1_sig)}

############################### BWrite writes out gene clusters from a BEAM-derived heatmap object in .csv format ######################
Bwrite <- function(beamhm,branchpoint){
  beamdf <- data.frame(cluster = beamhm$annotation_row$Cluster,gene = row.names(beamhm$annotation_row),stringsAsFactors = FALSE)
  for(x in unique(beamdf$cluster)){
    filename <- paste0('Basal_BEAM_',branchpoint,'_C',x,'.csv')
    write.csv(beamdf[which(beamdf$cluster==x),] file = filename)}}

########################################################################################################################################

## Do BEAM on all 5 branchpoints. 
BEAM_1 <- BEAMer(bmono,1)
BEAM_2 <- BEAMer(bmono,2)
BEAM_3 <- BEAMer(bmono,3)
BEAM_4 <- BEAMer(bmono,4)
BEAM_5 <- BEAMer(bmono,5)

## Plot the significantly branch dependent genes, adjust expression cluster numbers as necessary for each BEAM. 
# Return each BEAM plot as heatmaps for easy extraction of gene expression modules. 
beamgenes_1 <- plot_genes_branched_heatmap(bmono[row.names(BEAM_1),], branch_point = 1, cores = 6, use_gene_short_name = T, show_rownames = F, return_heatmap = T, num_clusters = 4)
beamgenes_2 <- plot_genes_branched_heatmap(bmono[row.names(BEAM_2),], branch_point = 2, cores = 6, use_gene_short_name = T, show_rownames = F, return_heatmap = T, num_clusters = 4)
beamgenes_3 <- plot_genes_branched_heatmap(bmono[row.names(BEAM_3),], branch_point = 3, cores = 6, use_gene_short_name = T, show_rownames = F, return_heatmap = T, num_clusters = 4)
beamgenes_4 <- plot_genes_branched_heatmap(bmono[row.names(BEAM_4),], branch_point = 4, cores = 6, use_gene_short_name = T, show_rownames = F, return_heatmap = T, num_clusters = 4)
beamgenes_5 <- plot_genes_branched_heatmap(bmono[row.names(BEAM_5),], branch_point = 5, cores = 6, use_gene_short_name = T, show_rownames = F, return_heatmap = T, num_clusters = 4)

## Write out gene clusters for each BEAM. 
Bwrite(beamgenes_1,1)
Bwrite(beamgenes_2,2)
Bwrite(beamgenes_3,3)
Bwrite(beamgenes_4,4)
Bwrite(beamgenes_5,5)
