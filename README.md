# Whelan_scRNA_Esophagus_May20
Seurat and Monocle commands for the Whelan Lab scRNA analysis on the murine esophageal epithelium.

This repository contains the R command lines for the analyses that were generated for the Whelan Lab paper <>. Raw data and matrices are available in GEO <>.
Analyses were done chiefly using Seurat(v3) and Monocle(v2). Additional libraries are specified throughout. 
There are 4 walkthroughs:
1. Initial_clustering_Seurat (Initial pre-processing and scRNA-Seq QC using Seurat, as well as some early visualizations of the data)
2. Monocle_PT_Epi (Pseudotime command lines using Monocle for all epithelial cells, along with BEAM)
3. Monocle_PT_Basal (Pseudotime for just basal cells, and BEAMs associated with it)
4. Downstream_Seurat (Additional Seurat analyses for data validation, miscellaneous differential expression tests, and marker validations)
