## Load Seurat and the epithelial file. ggplot2 and EnhancedVolcano will be used for some visualization, as well.
library(ggplot2)
library(EnhancedVolcano)
library(Seurat)
load('ovy.epi')

### AGING COMPARISONS

## Calculate proportional differences between clusters in aging. And find log odds. 
# Set the identities as 'NamedTypes', as these are the defined clusters.
Idents(ovy.epi) <- 'NamedTypes'
# Create a dataframe with proportions of each cluster in both age groups, and calculate log oddds.
cpt <- as.data.frame.matrix(prop.table(table(Idents(ovy.epi), ovy.epi$Age), margin = 2))
cpt$Odds <- cpt$Aged/cpt$Young
cpt$LogOdds <- log2(cpt$Odds)

## Find differentially expressed genes in each cluster between the age groups. 
# Calculate fold changes for ALL genes in each cluster.
B1ovy <- FindMarkers(ovy.epi, subset.ident = 'Basal-1',group.by = 'Age', ident.1 = 'Aged', logfc.threshold = 0, min.pct = 0, test.use = 'MAST')
B2ovy <- FindMarkers(ovy.epi, subset.ident = 'Basal-2',group.by = 'Age', ident.1 = 'Aged', logfc.threshold = 0, min.pct = 0, test.use = 'MAST')
B3ovy <- FindMarkers(ovy.epi, subset.ident = 'Basal-3',group.by = 'Age', ident.1 = 'Aged', logfc.threshold = 0, min.pct = 0, test.use = 'MAST')
B4ovy <- FindMarkers(ovy.epi, subset.ident = 'Basal-4',group.by = 'Age', ident.1 = 'Aged', logfc.threshold = 0, min.pct = 0, test.use = 'MAST')
B5ovy <- FindMarkers(ovy.epi, subset.ident = 'Basal-5',group.by = 'Age', ident.1 = 'Aged', logfc.threshold = 0, min.pct = 0, test.use = 'MAST')
Iovy <- FindMarkers(ovy.epi, subset.ident = 'Suprabasal',group.by = 'Age', ident.1 = 'Aged', logfc.threshold = 0, min.pct = 0, test.use = 'MAST')
D1ovy <- FindMarkers(ovy.epi, subset.ident = 'Superficial-1',group.by = 'Age', ident.1 = 'Aged', logfc.threshold = 0, min.pct = 0, test.use = 'MAST')
D2ovy <- FindMarkers(ovy.epi, subset.ident = 'Superficial-2',group.by = 'Age', ident.1 = 'Aged', logfc.threshold = 0, min.pct = 0, test.use = 'MAST')
# Create a master dataframe that combines all the fold changes in each cluster.
B1ovy$Cluster <- rep('Basal-1',length(rownames(B1ovy)))
B2ovy$Cluster <- rep('Basal-2',length(rownames(B2ovy)))
B3ovy$Cluster <- rep('Basal-3',length(rownames(B3ovy)))
B4ovy$Cluster <- rep('Basal-4',length(rownames(B4ovy)))
B5ovy$Cluster <- rep('Basal-5',length(rownames(B5ovy)))
Iovy$Cluster <- rep('Suprabasal',length(rownames(Iovy)))
D1ovy$Cluster <- rep('Superficial-1',length(rownames(D1ovy)))
D2ovy$Cluster <- rep('Superficial-2',length(rownames(D2ovy)))
# Assign "significant" or "non-significant" for each fold change given adjusted p-value.
jitdf <- as.data.frame(rbind(B1ovy,B2ovy,B3ovy,B4ovy,B5ovy,Iovy,D1ovy,D2ovy))
sig <- c()
for (x in jitdf$p_val_adj){if (x >= 0.05){sig <- append(sig,'Non-Significant')} else {sig <- append(sig,'Significant')}}
# Add significance classifications to the dataframe.
jitdf$Significance <- sig
# Relevel the clusters, just in case. 
jitdf$Cluster <- factor(jitdf$Cluster, levels = c('Basal-1','Basal-2','Basal-3','Basal-4','Basal-5','Intermediate','Differentiating-1','Differentiating-2'))
# Plot.
ggplot(jitdf, aes(x=Cluster, y=avg_logFC, color = Significance)) + geom_jitter() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab('Average Natural log FC') + scale_color_manual(values=c('blue','red'))

## Bulk differential expression test of old versus young.
OvY <- FindMarkers(ovy.epi, group.by = 'Age', ident.1 = 'Aged', test.use = 'MAST')
# Get just the significant DEGs.
OvY <- OvY[OvY$p_val_adj<0.05,]
# Write out DEGs.
write.csv(as.data.frame(OvY), file = 'OvY_DEGs.csv')
# Visualize DEGs' fold changes and p-values with a volcano plot.
EnhancedVolcano(OvY, rownames(OvY), 'avg_logFC','p_val_adj', FCcutoff = 0, pCutoff = 0.05, drawConnectors = F, transcriptPointSize = 3, transcriptLabSize = 3, xlab = 'Average natural log FC')

### EXPRESSION OF CANONICAL MARKERS

## Visualize canonical basal markers.
basalmarkers <- c('Mki67','Cd34','Nt5e', 'Itgb1','Itgb4','Itga6','Trp63','Sox2', 'Krt15','Krt14','Krt5')
DotPlot(ovy.epi, features = basalmarkers, group.by = 'NamedTypes', assay = 'RNA', cols = c('blue','red'))+ RotatedAxis()

## Visualize canonical superficial markers.
diffmarkers <- c('Ivl','Flg','Fst','Lor','Krt4','Krt13','Krtdap')
DotPlot(ovy.epi, features = diffmarkers, group.by = 'NamedTypes', assay = 'RNA', cols = c('blue','red'))+ RotatedAxis()

### FIND CLUSTER-EXCLUSIVE MARKERS

## These markers will be ranked by how lowly they're expressed in clusters other than the one in question.

# Run DEG tests on all clusters, but specify lower minimal detection.
bmarks <- FindAllMarkers(ovy.int, assay = 'RNA', test.use = 'MAST', min.pct = 0)
# Rank by lowest expression outside of it.
heatmarkers <-bmarks[which(bmarks$p_val_adj<0.05),] %>% group_by(cluster) %>% top_n(n = -5, wt = pct.2)
bavex <- AverageExpression(ovy.epi, return.seurat = T)
# Visualize.
DoHeatmap(object = bavex, features = heatmarkers$gene, assay = 'RNA', label = F, draw.lines = F)+scale_fill_gradientn(colors = c("blue", "white", "red"))

