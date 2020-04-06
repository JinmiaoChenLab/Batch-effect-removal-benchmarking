
# Author : Marion Chevrier
# Date : 01/10/2019
# Proj : scMerge pipeline

#' scMerge pipeline
#' 
#' @param expr_mat : expression matrix, genes by cells
#' @param metadata : dataframe with metadata
#' @param batch_label : labels in metadata for batch (default "batchlb")
#' @param celltype_label : labels in metadata for celltype (default "CellType")
#' @param filter_genes : logical for gene filtering (default FALSE)
#' @param filter_cells : logical for cell filtering (default FALSE)
#' @param min_cells : min cells for gene filtering (default 10)
#' @param min_genes : min genes for cell filtering (default 300)
#' @param cosineNorm : logical for cosine normalization (default FALSE)
#' @param genelist_ENS : logical for predefined vector of stably expressed genes (SEG) from scMerge package expressed as ensemblGeneID (otherwise expressed as gene name)
#' @param genelist_human : logical for predefined vector of stably expressed genes (SEG) from scMerge package expressed as human genes (otherwise expressed as mouse genes)
#' @param kmeansK : a vector indicates the kmeans's K for each batch
#' @param replicate_prop : a number indicating the ratio of cells that are included in pseudo-replicates, ranges from 0 to 1 (default 0.5)
#' @param npcs : number of principal components (default 20)
#' @param plotout_dir : path to saving plot files (default current directory)
#' @param saveout_dir : path to saving output files (default current directory)
#' @param save_obj : logical to save Seurat object (default TRUE)

scMerge_preprocess <- function(expr_mat, metadata, batch_label = "batchlb",
                               filter_genes = F, filter_cells = F,
                               min_cells = 10, min_genes = 300, cosineNorm= FALSE)
{
  
  ##########################################################
  # preprocessing
  
  if(filter_genes == F) {
    min_cells = 0
  }
  if(filter_cells == F) {
    min_genes = 0
  }
  
  num.genes <- colSums(expr_mat > 0)
  num.mol <- colSums(expr_mat)
  cells.use <- names(x = num.genes[which(x = num.genes > min_genes)])
  expr_mat <- expr_mat[, cells.use]
  
  genes.use <- rownames(expr_mat)
  num.cells <- rowSums(expr_mat > 0)
  genes.use <- names(x = num.cells[which(x = num.cells > min_cells)])
  expr_mat <- expr_mat[genes.use, ]
  
  metadata <- metadata[colnames(expr_mat),]
  
  ##########################################################
  # log / cosineNorm
  
  # Scale the data
  # note : as.matrix() does not work if matrix has more than 200000 rows  
  
  chunk <- function(x,n){
    vect <- c(1:x)
    num <- ceiling(x/n)
    split(vect,rep(1:num,each=n,len=x))
  }
  
  log_fun <- function(data){
    if(dim(data)[2]>200000){
      l <- lapply(chunk(ncol(data),200000), function(i){log2(data[,i]+1)})
      res <- do.call(cbind,l)
    } else {
      res <- log2(data+1)
    }
    return(res)
  }
  
  if(cosineNorm){
    cosdata_list <- lapply(levels(metadata[,batch_label]),function(l){
      mat <- expr_mat[,colnames(expr_mat) %in% rownames(metadata[metadata[,batch_label]==l,])]
      cosdata <- cosineNorm(as.matrix(log_fun(mat)))
      colnames(cosdata) <- colnames(mat)
      rownames(cosdata) <- rownames(mat)
      return(cosdata)
    })
    data_log <- do.call(cbind,cosdata_list)
    data_log <- data_log[,rownames(metadata)]
  } else {
    data_log=log_fun(expr_mat)
  }
  
  ##########################################################
  # SingleCellExperiment object
  
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(expr_mat), logcounts=as.matrix(data_log)), colData = metadata)

  return(sce)
}

call_scMerge <- function(sce, batch_label, celltype_label,
                         seg = NULL, 
                         kmeansK = kmeansK, marker = NULL, 
                         replicate_prop=0.5,npcs=20,
                         plotout_dir = "", saveout_dir = "", 
                         outfilename_prefix = "", 
                         visualize = T, save_obj = T)
{
  
  #plotting
  k_seed = 10
  tsne_perplex = 30
  tsneplot_filename = '_scMerge_tsne'
  obj_filename = "_scMerge_sobj"
  pca_filename = "_scMerge_pca"
  
  ##########################################################
  #run
  
  t1 = Sys.time()

  scMerge_res <- scMerge(
    sce_combine = sce, 
    ctl = seg,
    kmeansK = kmeansK,
    marker = marker,
    assay_name = "scMerge_res",
    replicate_prop = replicate_prop,
    cell_type = NULL, # unsupervised
    verbose=T,
    fast_svd=TRUE,
    parallel = FALSE)
  
  t2 = Sys.time()
  print(t2-t1)
  
  ##########################################################
  # post processing and save
  
  # Get PCA coordinates on normalized data
  normalized_data <- as.data.frame(scMerge_res@assays$data$scMerge_res) 
  normalized_data_t <- t(normalized_data)
  pca_mat <- stats::prcomp(normalized_data_t, rank = npcs, retx=TRUE, center = TRUE, scale. = FALSE)
  pca_mat <- as.data.frame(pca_mat$x)
  
  pca_mat[rownames(metadata), batch_label] <- metadata[, batch_label]
  pca_mat[rownames(metadata), celltype_label] <- metadata[ , celltype_label]
  write.table(pca_mat, file=paste0(saveout_dir, outfilename_prefix, pca_filename, ".txt"), quote=F, sep='\t', row.names = T, col.names = NA)
  
  if(save_obj) {
    saveRDS(scMerge_res, file=paste0(saveout_dir,outfilename_prefix,obj_filename,".RDS"))
  }
  
  if (visualize) {
    
    ##########################################################
    #preparing plots
    
    set.seed(k_seed)
    out_tsne <- Rtsne::Rtsne(pca_mat, perplexity = tsne_perplex, pca = TRUE, check_duplicates=FALSE)
    tsne_df <- as.data.frame(out_tsne$Y)
    
    ##########################################################
    #tSNE plot
    
    rownames(tsne_df) <- rownames(pca_mat)
    colnames(tsne_df) <- c('tSNE_1', 'tSNE_2')
    tsne_df[rownames(metadata), batch_label] <- metadata[, batch_label]
    tsne_df[rownames(metadata), celltype_label] <- metadata[ , celltype_label]
    
    p1 <- ggplot(tsne_df, aes_string(x = 'tSNE_1', y = 'tSNE_2', color = batch_label)) + 
      geom_point(alpha = 0.6) + theme_bw() + ggtitle('scMerge') + 
      theme(legend.title = element_text(size=17), 
            legend.key.size = unit(1.1, "cm"),
            legend.key.width = unit(0.5,"cm"), 
            legend.text = element_text(size=14), 
            plot.title = element_text(color="black", size=20, hjust = 0.5))
    
    p2 <- ggplot(tsne_df, aes_string(x = 'tSNE_1', y = 'tSNE_2', color = celltype_label)) + 
      geom_point(alpha = 0.6) + theme_bw() + ggtitle('scMerge') + 
      theme(legend.title = element_text(size=17), 
            legend.key.size = unit(1.1, "cm"),
            legend.key.width = unit(0.5,"cm"), 
            legend.text = element_text(size=14), 
            plot.title = element_text(color="black", size=20, hjust = 0.5))
    
    png(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p1, p2))
    dev.off()
    
    pdf(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".pdf"),width=15,height=7,paper='special')
    print(plot_grid(p1, p2))
    dev.off()
  }
  
  return(scMerge_res)
}

