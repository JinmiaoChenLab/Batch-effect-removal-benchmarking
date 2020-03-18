
# Author : Kok Siong Ang 
# Date : 18/09/2019
# Proj : MNN Correct pipeline

#' MNN Correct pipeline
#' 
#' @param expr_mat_mnn : expression matrix for seurat object, genes by cells
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
#' @param npcs : number of principal components (default 20)
#' @param batch_label : labels in metadata for batch (default "batchlb")
#' @param celltype_label : labels in metadata for celltype (default "CellType")
#' 
#' @param plotout_dir : path to saving plot files (default current directory)
#' @param saveout_dir : path to saving output files (default current directory)
#' @param save_obj : logical to save Seurat object (default TRUE)

mnn_preprocess <- function(expr_mat_mnn, metadata, 
                              filter_genes = T, filter_cells = T,
                              normData = T, Datascaling = T, regressUMI = F, 
                              min_cells = 10, min_genes = 300, norm_method = "LogNormalize", scale_factor = 10000, 
                              b_x_low_cutoff = 0.0125, b_x_high_cutoff = 3, b_y_cutoff = 0.5, 
                              numVG = 300, npcs = 20, 
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

  b_seurat <- CreateSeuratObject(raw.data = expr_mat_mnn, meta.data = metadata, project = "MNNcorrect_benchmark", 
                                 min.cells = min_cells, min.genes = min_genes)
  if (normData) {
    b_seurat <- NormalizeData(object = b_seurat, normalization.method = norm_method, scale.factor = scale_factor)
  } else {
    b_seurat@data = b_seurat@raw.data
  }

  if(regressUMI && Datascaling) {
    b_seurat <- ScaleData(object = b_seurat, vars.to.regress = c("nUMI"))  # in case of read count data
  } else if (Datascaling) { # default option
    b_seurat <- ScaleData(object = b_seurat)
  } else {
    b_seurat@scale.data <- b_seurat@data
  }

  hvg_all <- rownames(b_seurat@scale.data)
  
  #split data into two batches
  meta_data <- b_seurat@meta.data

  b1_selected <- row.names(meta_data)[which (meta_data$batch==1)]
  b2_selected <- row.names(meta_data)[which (meta_data$batch==2)]
  b1_meta <- meta_data[which (meta_data$batch==1),]
  b2_meta <- meta_data[which (meta_data$batch==2),]

  data_filtered <- b_seurat@scale.data[hvg_all,]
  dim(data_filtered)
  b1_data <- data_filtered[,b1_selected]
  b2_data <- data_filtered[,b2_selected]
  b1_meta <- b1_meta[b1_selected,]
  b2_meta <- b2_meta[b2_selected,]

  return(list( b1_data, b2_data, b1_meta, b2_meta ))
}

call_mnn <- function(batch_data, batch_label, celltype_label, npcs = 20, 
                         plotout_dir = "", saveout_dir = "", 
                         outfilename_prefix = "", 
                         visualize = T, save_obj = T)
{

  #plotting
  k_seed = 10
  tsne_perplex = 30
  tsneplot_filename = "_mnncorrect_tsne"
  mnn_out_filename = "_mnncorrect_out"
  metadata_out_filename = "_mnncorrect_metadata_out"
  pca_filename = "_mnncorrect_pca"

  b1_data <- batch_data[[1]]
  b2_data <- batch_data[[2]]
  b1_meta <- batch_data[[3]]
  b2_meta <- batch_data[[4]]

  metadata <- rbind(b1_meta, b2_meta)

  ##########################################################
  #run

  t1 = Sys.time()

  out_mnn_total <- mnnCorrect(as.matrix(b1_data), as.matrix(b2_data), 
                              k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE)

  t2 = Sys.time()
  print(t2-t1)

  ##########################################################
  #save

  saveRDS(out_mnn_total, file=paste0(saveout_dir, outfilename_prefix, mnn_out_filename, ".RDS"))

  mnn_cor_t <- as.matrix(t(do.call(cbind, out_mnn_total$corrected)))
  rownames(mnn_cor_t) <- colnames(cbind(b1_data,b2_data))
  cells_use <- rownames(mnn_cor_t)
  metadata <- metadata[cells_use,]
  #write.table(mnn_cor_t, file=paste0(saveout_dir, outfilename_prefix, mnn_out_filename, "_t.txt"), quote=F, sep="\t", row.name=TRUE, col.name=NA)

  saveRDS(metadata, file = paste0(saveout_dir, outfilename_prefix, metadata_out_filename, ".RDS"))
  #write.table(metadata, file=paste0(saveout_dir, outfilename_prefix, metadata_out_filename, ".txt"), quote=F, sep="\t", row.name=TRUE, col.name=NA)

  pca_mat <- prcomp(mnn_cor_t, rank=npcs)
  pca_mat <- as.data.frame(pca_mat$x)

  pca_mat[rownames(metadata), batch_label] <- metadata[ , batch_label]
  #pca_mat[rownames(metadata), 'batch'] <- metadata[ , "batch"]
  pca_mat[rownames(metadata), celltype_label] <- metadata[ , celltype_label]
  write.table(pca_mat, file=paste0(saveout_dir, outfilename_prefix, pca_filename, ".txt"), quote=F, sep='\t', row.names = T, col.names = NA)

  if (visualize) {

    set.seed(10)
    #out_tsne <- Rtsne(mnn_cor_t, check_duplicates=FALSE, perplexity = 30)
    out_tsne <- Rtsne(pca_mat, check_duplicates=FALSE, perplexity = 30)
    tsne_df<- as.data.frame(out_tsne$Y)

    rownames(tsne_df) <- colnames(cbind(b1_data,b2_data))
    dim(tsne_df)
    colnames(tsne_df) <- c('tSNE_1', 'tSNE_2')
    tsne_df$batchlb <- metadata[ , batch_label]
    #tsne_df$batch <- metadata[ , "batch"]
    tsne_df$CellType <- metadata[ , celltype_label]

    p01 <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2, colour = batchlb)) + geom_point(alpha = 0.6) + theme_bw()
    p01 <- p01 + labs(x='tSNE_1',y='tSNE_2',title='MNNcorrect')
    p01 <- p01 + theme(legend.title = element_text(size=17), 
                       legend.key.size = unit(1.1, "cm"),
                       legend.key.width = unit(0.5,"cm"), 
                       legend.text = element_text(size=14), 
                       plot.title = element_text(color="black", size=20, hjust = 0.5))

    p02 <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2, colour = CellType)) + geom_point(alpha = 0.6) + theme_bw()
    p02 <- p02 + labs(x='tSNE_1',y='tSNE_2',title='MNNcorrect')
    p02 <- p02 + theme(legend.title = element_text(size=17), 
                       legend.key.size = unit(1.1, "cm"),
                       legend.key.width = unit(0.5,"cm"), 
                       legend.text = element_text(size=14), 
                       plot.title = element_text(color="black", size=20, hjust = 0.5))

    png(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p01, p02))
    dev.off()

    pdf(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".pdf"),width=15,height=7,paper='special')
    print(plot_grid(p01, p02))
    dev.off()
  }

}







