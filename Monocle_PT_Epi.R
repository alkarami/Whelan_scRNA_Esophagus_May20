## Packages.
library(Seurat)
library(monocle)
library(cowplot)

############## newimport imports Seurat v3 object as Monocle2. Credits to Lina Kroehling (https://www.biostars.org/u/52871/). ##########
newimport <- function(otherCDS, import_all = T) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- otherCDS@assays$RNA@counts
    
    if(class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }
    
    pd <- tryCatch( {
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, 
    #warning = function(w) { },
    error = function(e) { 
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      
      message("This Seurat object doesn't provide any meta data");
      pd
    })
    
    # remove filtered cells from Seurat
    if(length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]  
    }
    
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- 0
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    valid_data <- data[, row.names(pd)]
    
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    
    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
        
      } else {
        # mist_list <- list(ident = ident) 
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    }
    
    if(1==1) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@assays$RNA@var.features)
      
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
    
  } else if (class(otherCDS)[1] == 'SCESet') {
    requireNamespace("scater")
    
    message('Converting the exprs data in log scale back to original scale ...')    
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
    
    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if("is.expr" %in% slotNames(otherCDS))
      lowerDetectionLimit <- otherCDS@is.expr
    else 
      lowerDetectionLimit <- 1
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    if(import_all) {
      # mist_list <- list(iotherCDS@sc3,
      #                   otherCDS@reducedDimension)
      mist_list <- otherCDS 
      
    } else {
      mist_list <- list()
    }
    
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    # monocle_cds@auxClusteringData$sc3 <- otherCDS@sc3
    # monocle_cds@auxOrderingData$scran <- mist_list
    
    monocle_cds@auxOrderingData$scran <- mist_list
    
  } else {
    stop('the object type you want to export to is not supported yet')
  }
  
  return(monocle_cds)
}

################################## monocon function is a wrapper to convert and process Seurat data ####################################
monocon <- function(o){
  # Convert OvY into Monocle2 with newimport
  all.mono <- newimport(o, import_all = T)
  # Pre-process
  all.mono <- estimateSizeFactors(all.mono)
  all.mono <- estimateDispersions(all.mono)
  ## FOLLOW DPFEATURE PROCEDURE
  
  # Preprocessing
  all.mono <- detectGenes(all.mono, min_expr = 0.1)
  fData(all.mono)$use_for_ordering <-
    fData(all.mono)$num_cells_expressed > 0.05 * ncol(all.mono)
  
  # Run tSNE, perform dimensionality reduction and regress out batch effects from samples.
  all.mono <- reduceDimension(all.mono,
                              norm_method = 'log',
                              reduction_method = 'tSNE',
                              residualModelFormulaStr = '~orig.ident',
                              verbose = T)
                              
  # Get DEGs that distinguish each cluster
  all.mono <- clusterCells(all.mono, verbose = F)
  clustering_DEG_genes <-
    differentialGeneTest(all.mono,
                         fullModelFormulaStr = '~Cluster', 
                         cores = 6)
  
  # Use top 2k genes to order trajectory
  all.mono_ordering_genes <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:2000]
  
  all.mono <-
    setOrderingFilter(all.mono,
                      ordering_genes = all.mono_ordering_genes)
  
  all.mono <-
    reduceDimension(all.mono, reduction_method = "DDRTree")
  
  all.mono <-
    orderCells(all.mono)
  return(all.mono)
}

## Load the Seurat object 'ovy.epi', the integrated epithelial dataset. Change assay to integrated.
load('ovy.epi')
DefaultAssay(ovy.epi) <- 'integrated'

## Convert ovy.epi into a Monocle object for pseudotime analyses.
wtmono <- monocon(ovy.epi)

## Visualize the pseudotime trajectory and the localization of Mki67, a proliferation marker. 
p1 <- plot_cell_trajectory(wtmono, color_by = 'NamedTypes')
p2 <- plot_cell_trajectory(wtmono, markers = 'Mki67', use_color_gradient = T)
plot_grid(p1,p2)

## Root the pseudotime at the point of proliferation (State 1).
wtmono <- orderCells(wtmono, root_state = 1)

## Identify significantly branch-dependent differentially expressed genes at the branchpoint with BEAM.
BEAM_1 <- BEAM(wtmono, branch_point = 1, cores = 4)
BEAM_1 <- BEAM_1[order(BEAM_1$qval),]
BEAM_1 <- BEAM_1[,c("gene_short_name", "pval", "qval")]
BEAM_1_sig <- BEAM_1[which(BEAM_1$qval<0.05),]

## Plot the genes, and cluster the expression modules. Save the gene clusters for pathway analyses.
beamgenes <- plot_genes_branched_heatmap(wtmono[row.names(BEAM_1_sig),], branch_labels = c("Differentiation", "Quiescence"),
                                         branch_point = 1, cores = 4, use_gene_short_name = T, show_rownames = F, return_heatmap = T, num_clusters = 7)

beamdf <- data.frame(cluster = beamgenes$annotation_row$Cluster,genes = row.names(beamgenes$annotation_row),stringsAsFactors = FALSE)
for(x in unique(beamdf$cluster)){
  filename <- paste0('AllBeam_C',x,'.csv')
  write.csv(beamdf[which(beamdf$cluster==x),] file = filename)
}

## Save the object.
save(wtmono, file = 'wtmono')
