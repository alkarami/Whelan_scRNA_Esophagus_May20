## Packages.
library(Seurat)

## Get inDrop files. 
all_files <- list.files('WT_Matrices/.')

## Initialize file list and upload all files. 
seurat_list <- c()
for (x in all_files){
  y <- paste0('WT_Matrices/',x)
  z <- read.table(file = y, sep = '\t', header = T, row.names = 1, as.is = T)
  seurat_list <- append(seurat_list,CreateSeuratObject(z, assay = "RNA",
                                                       min.cells = 0, min.features = 0, names.field = 1,
                                                       names.delim = "_", meta.data = NULL, project = x))
}

## QC and pre-process for all files.
# Pre-process: Add metadata to each file, normalize the matrices with log1p transformation (NormalizeData()), calculate percentage
# of transcripts expressing mitochondrial genes as percent of total cell expression (PercentageFeatureSet()), subset cells with
# less than 250 or more than 2500 unique transcripts and finally find 2000 most variable features in the matrix. For a permissive 
# initial clustering and our interest in mitochondrial function, high mitochondria-expressing cells are not yet subsetted out.

aged <- merge(seurat_list[[1]], unlist(seurat_list[2:6]), merge.data = T, project = 'Aged')
age <- sample(c("Aged"),size = length(colnames(aged)),replace = TRUE)
names(age) <- colnames(aged)
aged <- AddMetaData(aged,age,col.name = 'Age')
aged <- subset(aged, subset = nFeature_RNA > 250 & nFeature_RNA < 2500)
aged <- NormalizeData(aged)
aged <- FindVariableFeatures(aged, selection.method = "vst", nfeatures = 2000)

young <- merge(seurat_list[[7]], unlist(seurat_list[8:12]), merge.data = T, project = 'Young')
age <- sample(c("Young"),size = length(colnames(young)),replace = TRUE)
names(age) <- colnames(young)
young <- AddMetaData(young,age,col.name = 'Age')
young <- subset(young, subset = nFeature_RNA > 250 & nFeature_RNA < 2500)
young <- NormalizeData(young)
young <- FindVariableFeatures(young, selection.method = "vst", nfeatures = 2000)

## Get anchors & integrate data according to the anchors.
int.anchors <- FindIntegrationAnchors(object.list = list(aged,young), dims = 1:30)
ovy.int <- IntegrateData(anchorset = int.anchors, dims = 1:30)
DefaultAssay(ovy.int) <- "integrated"
## Scale the integrated data matrix.
ovy.int <- ScaleData(ovy.int, verbose = FALSE)
## Perform PCA and produce the first 20 principal components.
ovy.int <- RunPCA(ovy.int, npcs = 30, verbose = FALSE)
## Perform UMAP from the resulting PCA PC's.
ovy.int <- RunUMAP(ovy.int, reduction = "pca", dims = 1:20)
## Construct SNN to calculate clusters of cells. 
ovy.int <- FindNeighbors(ovy.int, dims = 1:20)
ovy.int <- FindClusters(ovy.int, resolution = 0.5)

## Visualize cells in UMAP.
DimPlot(ovy.int, split.by = 'Age')

## Perform post-clustering QC.
# Unique transcripts.
VlnPlot(ovy.int, features = c('nFeature_RNA'))
# Mitochondrial expression. Most of the distribution seems to stay around the 0-20% mark. Subset out cells with >20% percentage.
MTgenes <- c('Atp8','Atp6','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4L','ND4','ND5','ND6','Rnr2')
MTgenes <- MTgenes[which(MTgenes %in% rownames(ovy.int))]
ovy.int[['P.Mito']] <- PercentageFeatureSet(ovy.int, features = MTgenes)
VlnPlot(ovy.int, features = c('P.Mito'))
ovy.int <- subset(ovy.int, subset = P.Mito < 20)
# Sample-wise cell distribution in clusters. 
DimPlot(ovy.int, group.by = 'orig.ident')

## Cluster 7 is likely dead or dying cells, given uniquely elevated mitochondrial content and low feature expression. Subset out.
ovy.int <- subset(ovy.int, subset = c(0,1,2,3,4,5,6,8,9,10))

## Save clean clusters.
save(ovy.int, file = 'ovy.int')

## Identify strong UMAP contributor genes.
DimHeatmap(ovy.int, reduction = 'umap', projected = T, dims = 1, cells = 1000, balanced = T, fast = F)+scale_fill_gradientn(colors = c("blue", "white", "red"))

## Change Seurat assay into un-integrated data. 
DefaultAssay(ovy.int) <- 'RNA'

## Visualize expression of basal and superficial markers that were identified as strong UMAP contributors.
FeaturePlot(ovy.int, features = c('Krt5','Krtdap'))

## Rename clusters based on marker expression.
ovy.int$NamedTypes <- revalue(ovy.epi$seurat_clusters, replace = c('0' = 'Basal-1','1' = 'Suprabasal','3' = 'Basal-2','4' = 'Superficial-2','5' = 'Basal-3','6' = 'Basal-4','8' = 'Basal-5','9' = 'Immune','10' = 'Fibroblast'))

## Identify cluster markers, and visualize the top 5. Write out the differentially expressed genes into file for IPA analyses. 
bmarks <- FindAllMarkers(ovy.int, assay = 'RNA', test.use = 'MAST')
write.csv(as.data.frame(bmarks), file = 'ClusterWise_DEGs_Whelan_scEso.csv')
heatmarkers <-bmarks[which(bmarks$p_val_adj<0.05),] %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
bavex <- AverageExpression(ovy.int, return.seurat = T)
DoHeatmap(object = bavex, features = heatmarkers$gene, assay = 'RNA', label = F, draw.lines = F, group.colors = c('Basal-1' = 'steelblue', 'Basal-2' = 'deeppink2', 'Basal-3' = 'green3', 'Basal-4' = 'orange', 'Basal-5' = 'red', 'Intermediate' = 'black','Differentiating-1' = 'darkviolet','Differentiating-2'='gold1','Fibroblast' = 'grey55','Immune'='grey') )+scale_fill_gradientn(colors = c("blue", "white", "red"))

## Assess the cell cycle phase of each cell using the expression of cell cycle related genes. 
# Define genes that contribute to the S phase.
sgenes <- c()
for (x in cc.genes$s.genes){
  sgenes <- append(sgenes,paste(substring(x,1,1),tolower(substring(x,2)), sep = ''))
}
# Define genes that contribute to the G2/M phases.
g2mgenes <- c()
for (x in cc.genes$g2m.genes){
  g2mgenes <- append(g2mgenes,paste(substring(x,1,1),tolower(substring(x,2)), sep = ''))
}
# Calculate S phase and G2/M scores for each cell. 
ovy.int <- CellCycleScoring(ovy.int, s.features = sgenes, g2m.features = g2mgenes, set.ident = F)
# Visualize cell cycle phases.
DimPlot(ovy.int, group.by = 'Phase')

## Subset out immune cells for epithelium-exclusive downstream analyses and save.
ovy.epi <- subset(ovy.int, idents = c('Basal-1','Basal-2','Basal-3','Basal-4','Basal-5','Suprabasal','Superficial-1','Superficial-2'))
save(ovy.epi, file = 'ovy.epi')
