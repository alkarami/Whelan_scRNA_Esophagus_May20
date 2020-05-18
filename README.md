# Whelan_scRNA_Esophagus_May20
This repository contains the bash command lines for processing the raw data and R command lines for the analyses that were generated for the Whelan Lab paper <>. Raw data and matrices are available in GEO <>.
Analyses were done chiefly using Seurat(v3) and Monocle(v2). Additional libraries are specified throughout. 

**There are 5 walkthroughs:**
1. PreProcess_and_Alignment (Bash lines for read deconvolution, genome alignment, feature counting, and matrix generation for downstream analyses)
2. Initial_clustering_Seurat (Initial pre-processing and scRNA-Seq QC using Seurat, as well as some early visualizations of the data)
3. Monocle_PT_Epi (Pseudotime command lines using Monocle for all epithelial cells, along with BEAM)

*Note that functions defined in Monocle_PT_Epi are used in Monocle_PT_Basal*

4. Monocle_PT_Basal (Pseudotime for just basal cells, and BEAMs associated with it)
5. Downstream_Seurat (Additional Seurat analyses for data validation, miscellaneous differential expression tests, and marker validations)

**Major requirements:** 

Bash: 

UMI-Tools (1.0.0): https://umi-tools.readthedocs.io/en/latest/

STAR (2.7): https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

Subread (2.0): http://subread.sourceforge.net/

R: 

Seurat (3.1.4): https://satijalab.org/seurat/

Monocle (2.16.0): http://cole-trapnell-lab.github.io/monocle-release/docs_mobile/
