
# Author : Kok Siong Ang 
# Date : 03/09/2019
# Proj : Seurat 2 pipeline

#' Seurat 2 pipeline
#' 
#' @param expr_mat_seurat : expression matrix for seurat object, genes by cells
#' @param metadata : dataframe with metadata for seurat object
#' @param filter_genes : logical for gene filtering (default TRUE)
#' @param filter_cells : logical for cell filtering (default TRUE)
#' @param normData : logical for expression matrix normalization (default TRUE)
#' @param Datascaling : logical for expression matrix scaling (default TRUE)
#' @param regressUMI : logical for UMI (default FALSE)
#' @param min_cells : min cells for gene filtering (default 10)
#' @param min_genes : min genes for cell filtering (default 300)
#' @param b_x_low_cutoff : Seurat FindVariableGenes, bottom cutoff on x-axis for identifying variable genes (default 0.0125)
#' @param b_x_high_cutoff : Seurat FindVariableGenes, top cutoff on x-axis for identifying variable genes (default 3)
#' @param b_y_cutoff : Seurat FindVariableGenes, top cutoff on y-axis for identifying variable genes (default 0.5)
#' @param scale_factor :  scale factor for normalization (default 10000)
#' @param numVG : number of variable genes for seurat (default 300)
#' @param nhvg : number of highly variable genes (default 2000)
#' @param batch_label : labels in metadata for batch (default "batchlb")
#' @param celltype_label : labels in metadata for celltype (default "CellType")
#' 
#' @param plotout_dir : path to saving plot files (default current directory)
#' @param saveout_dir : path to saving output files (default current directory)
#' @param save_obj : logical to save Seurat object (default TRUE)

seurat2_preprocess <- function(expr_mat_seurat, metadata, 
                              filter_genes = T, filter_cells = T,
                              normData = T, Datascaling = T, regressUMI = F, 
                              min_cells = 10, min_genes = 300, norm_method = "LogNormalize", scale_factor = 10000, 
                              b_x_low_cutoff = 0.0125, b_x_high_cutoff = 3, b_y_cutoff = 0.5, 
                              numVG = 300, nhvg = 2000, 
                              batch_label = "batchlb", celltype_label = "CellType")
{

  ##########################################################
  # preprocessing

  if(filter_genes == F) {
    min_cells = 0
  }
  if(filter_cells == F) {
    min_genes = 0
  }

  count = 0
  batch_list = list()
  hvg_union = list()
  for (label in unique(metadata[, batch_label])) {
    count = count +1
    selected <- metadata[, batch_label] == label
    metadata_batch <- metadata[selected ,]
    expr_mat_batch <- expr_mat_seurat[, selected]

    b_seurat <- CreateSeuratObject(raw.data = expr_mat_batch, meta.data = metadata_batch, 
                                   project = paste0("seurat2_",label), 
                                   min.cells = min_cells, min.genes = min_genes)

    if (normData) {
      b_seurat <- NormalizeData(object = b_seurat, normalization.method = "LogNormalize", scale.factor = scale_factor)
    } else {
      b_seurat@data = b_seurat@raw.data
    }

    b_seurat <- FindVariableGenes(object = b_seurat, do.plot = TRUE, mean.function = ExpMean, dispersion.function = LogVMR, 
                                  x.low.cutoff = b_x_low_cutoff, x.high.cutoff = b_x_high_cutoff, y.cutoff = b_y_cutoff)
    print(paste("Num var genes: ", length(x = b_seurat@var.genes)))
    if(length(x = b_seurat@var.genes) < numVG) {
      print(paste0("The number of highly variable genes in batch ", label, " is small, can change the threshold to gain more HVG genes"))
    }

    if(regressUMI && Datascaling) {
      b_seurat <- ScaleData(object = b_seurat, vars.to.regress = c("nUMI"))  # in case of read count data
    } else if (Datascaling) { # default option
      b_seurat <- ScaleData(object = b_seurat)
    }
    batch_list = c(batch_list, b_seurat)

    # Get the list of highly variable gens in 2 batches
    hvg_b <- rownames(x = head(x = b_seurat@hvg.info, n = nhvg))
    hvg_union <- union(x = hvg_union, y = hvg_b)
  }

  for (b_seurat in batch_list) {
    g <- rownames(b_seurat@scale.data)
    hvg_union <- intersect(hvg_union, g)
  }
  length(hvg_union)

  return(list(batch_list, hvg_union))
}

call_seurat2 <- function(batch_list, hvg_union, batch_label, celltype_label, npcs = 20, 
                         plotout_dir = "", saveout_dir = "", 
                         outfilename_prefix = "", 
                         visualize = T, save_obj = T)
{

  #plotting
  k_seed = 10
  tsne_perplex = 30
  tsneplot_filename = '_seurat2_tsne'
  umapplot_filename = '_seurat2_umap'
  obj_filename = "_seurat2_sobj"
  cca_filename = "_seurat2_cca"
  pca_filename = "_seurat2_pca"

  ##########################################################
  #run

  t1 = Sys.time()

  b_seurat = NULL
  if (length(batch_list) == 2) {
    b_seurat <- RunCCA(object = batch_list[[1]], object2 = batch_list[[2]], genes.use = hvg_union)
  } else {
    b_seurat = RunMultiCCA(batch_list, genes.use=hvg_union, num.ccs = npcs)
  }
  b_seurat <- AlignSubspace(object = b_seurat, reduction.type = "cca", grouping.var = batch_label, dims.align = 1:npcs)

  t2 = Sys.time()
  print(t2-t1)

  ##########################################################
  #post processing and save

  seurat2_res <- as.data.frame(b_seurat@dr$cca.aligned@cell.embeddings)
  seurat2_res$batchlb <- b_seurat@meta.data[, batch_label]
#  seurat2_res$batch <- b_seurat@meta.data[, ]
  seurat2_res$celltype <- b_seurat@meta.data[, celltype_label]
  write.table(seurat2_res, file=paste0(saveout_dir,outfilename_prefix,cca_filename,".txt"), quote=F, sep='\t', row.names = T, col.names = NA)

  mnn_corr_proc <- b_seurat@dr$cca.aligned@cell.embeddings[,1:npcs]
  tmp <- as.data.frame(mnn_corr_proc)
  #str(tmp)
  pca_mnn_corr <- prcomp( tmp , rank=npcs)

  pca_mnn_corr <- as.data.frame(pca_mnn_corr$x)
  pca_mnn_corr$batchlb <- b_seurat@meta.data[, batch_label]
#  pca_mnn_corr$batch <- b_seurat@meta.data[, ]
  pca_mnn_corr$celltype <- b_seurat@meta.data[, celltype_label]
  write.table(pca_mnn_corr, file=paste0(saveout_dir,outfilename_prefix,pca_filename,".txt"), quote=F, sep='\t', row.names = T, col.names = NA)

  if(save_obj) {
    saveRDS(b_seurat, file=paste0(saveout_dir,outfilename_prefix,obj_filename,".RDS"))
  }

  if (visualize) {

    ##########################################################
    #preparing plots

    b_seurat <- RunTSNE(b_seurat, reduction.use = "cca.aligned", dims.use = 1:npcs, k.seed = k_seed, do.fast = T, check_duplicates = FALSE, perplexity = tsne_perplex)
    b_seurat <- RunUMAP(b_seurat, reduction.use = "cca.aligned", dims.use = 1:npcs, k.seed = k_seed, do.fast = T)

    ##########################################################
    #tSNE plot

    p11 <- TSNEPlot(b_seurat, do.return = T, pt.size = 0.5, group.by = batch_label)
    p12 <- TSNEPlot(b_seurat, do.return = T, pt.size = 0.5, group.by = celltype_label)

    png(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p11, p12))
    dev.off()
  
    pdf(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".pdf"),width=15,height=7,paper='special') 
    print(plot_grid(p11, p12))
    dev.off()

    ##########################################################
    #UMAP plot

    p21 <- DimPlot(object = b_seurat, reduction.use = 'umap', group.by = batch_label)
    p22 <- DimPlot(object = b_seurat, reduction.use = 'umap', group.by = celltype_label)

    png(paste0(plotout_dir,outfilename_prefix,umapplot_filename,".png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p21, p22))
    dev.off()

    pdf(paste0(plotout_dir,outfilename_prefix,umapplot_filename,".pdf"),width=15,height=7,paper='special') 
    print(plot_grid(p21, p22))
    dev.off()

  }

  return(b_seurat)
}







