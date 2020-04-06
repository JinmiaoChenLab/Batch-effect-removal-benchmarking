
# Author : Marion Chevrier
# Date : 30/09/2019
# Proj : LIGER pipeline

#' LIGER pipeline
#' 
#' @param expr_mat : expression matrix, genes by cells
#' @param metadata : dataframe with metadata
#' @param batch_label : labels in metadata for batch (default "batchlb")
#' @param celltype_label : labels in metadata for celltype (default "CellType")
#' @param var.thresh : variance threshold used to identify variable genes. Genes with expression variance greater than threshold (relative to mean) are selected. (default 0.1)
#' @param k : number of factors (default 20)
#' @param nrep : number of restarts to perform factorization (default 3)
#' @param plotout_dir : path to saving plot files (default current directory)
#' @param saveout_dir : path to saving output files (default current directory)
#' @param save_obj : logical to save Seurat object (default TRUE)

liger_preprocess <- function(expr_mat, metadata, batch_label = "batchlb", var.thresh=0.1)
{
  
  data <- lapply(levels(metadata[,batch_label]), function(l){
    cells <- metadata$cell[metadata[,batch_label]==l]
    return(expr_mat[,cells])
  })
  names(data) <- levels(metadata[,batch_label])
  
  ##########################################################
  # preprocessing

  # Create LIGER object
  ligerex = createLiger(data) 
 
  # Normalize the data to account for different numbers of UMIs per cell 
  ligerex = normalize(ligerex)  
  
  # Select variable genes on each of the datasets separately, then takes the union
  ligerex = selectGenes(ligerex, var.thresh = var.thresh)
  
  # Scale the data
  # note : as.matrix() does not work if matrix has more than 200000 rows  
  
  dim_data <- sapply(data,function(l){dim(l)[1]})
  
  chunk <- function(x,n){
    vect <- c(1:x)
    num <- ceiling(x/n)
    split(vect,rep(1:num,each=n,len=x))
  }
  
  if(any(dim_data>200000)){
    ligerex@scale.data <- lapply(1:length(ligerex@norm.data),function(i) {liger:::scaleNotCenterFast(t(ligerex@norm.data[[i]][ligerex@var.genes,]))})
    ligerex@scale.data <- lapply(ligerex@scale.data,function(l){
      if(dim(l)[1]>200000){
        l2 <- lapply(chunk(nrow(l),200000), function(i){as.matrix(l[i,])})
        res <- do.call(rbind,l2)
      } else {
        res <- as.matrix(l)
      }
      return(res)
    })
    names(ligerex@scale.data) <- names(ligerex@norm.data)
    for (i in 1:length(ligerex@scale.data)) {
      ligerex@scale.data[[i]][is.na(ligerex@scale.data[[i]])] <- 0
      rownames(ligerex@scale.data[[i]]) <- colnames(ligerex@raw.data[[i]])
      colnames(ligerex@scale.data[[i]]) <- ligerex@var.genes
    }
    ligerex <- liger:::removeMissingObs(ligerex, slot.use = "scale.data", use.cols = F)
  } else {
    ligerex = scaleNotCenter(ligerex)
  }
  
  return(ligerex)
}


call_liger <- function(ligerex, metadata, batch_label, celltype_label,
                       k = 20, nrep = 3,
                       plotout_dir = "", saveout_dir = "", 
                       outfilename_prefix = "", 
                       visualize = T, save_obj = T)
{
  
  #plotting
  k_seed = 10
  tsne_perplex = 30
  tsneplot_filename = '_liger_tsne'
  obj_filename = "_liger_sobj"
  pca_filename = "_liger_pca"
  
  ##########################################################
  #run
  
  t1 = Sys.time()
  
  # perform the factorization using an alternating least squares algorithm
  ligerex = optimizeALS(ligerex, k = k, lambda = 5, nrep = nrep)
  # SNF clustering and quantile alignment
  ligerex = quantileAlignSNF(ligerex)
  
  t2 = Sys.time()
  print(t2-t1)
  
  ##########################################################
  #post processing and save
  
  ligerex_res <- as.data.frame(ligerex@H.norm)
  cells_use <- rownames(ligerex_res)
  metadata <- metadata[cells_use,]
  ligerex_res$batchlb <- metadata[, batch_label]
  ligerex_res$celltype <- metadata[, celltype_label]
  write.table(ligerex_res, file=paste0(saveout_dir,outfilename_prefix,pca_filename,".txt"), quote=F, sep='\t', row.names = T, col.names = NA)
  
  if(save_obj) {
    saveRDS(ligerex, file=paste0(saveout_dir,outfilename_prefix,obj_filename,".RDS"))
  }
  
  if (visualize) {
    
    ##########################################################
    #preparing plots
    
    ligerex = liger::runTSNE(ligerex, perplexity = tsne_perplex, method = "Rtsne", rand.seed = k_seed)
    
    ##########################################################
    #tSNE plot
    
    # Get coordinates 
    df <- data.frame(ligerex@tsne.coords)
    cells_use <- rownames(df)
    metadata <- metadata[cells_use,]
    
    ligerex@clusters <- as.factor(metadata[,batch_label])
    names(ligerex@clusters) <- rownames(metadata)
    p11 <- plotByDatasetAndCluster(ligerex, pt.size = 1, text.size = 0)
    
    ligerex@clusters <- as.factor(metadata[,celltype_label])
    names(ligerex@clusters) <- rownames(metadata)
    p12 <- plotByDatasetAndCluster(ligerex, pt.size = 1, text.size = 0)
    
    png(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p11, p12))
    dev.off()
    
    pdf(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".pdf"),width=15,height=7,paper='special') 
    print(plot_grid(p11, p12))
    dev.off()
    
  }
  
  return(ligerex)
}

