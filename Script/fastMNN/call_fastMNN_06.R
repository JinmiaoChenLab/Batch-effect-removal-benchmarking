
# Author : Kok Siong Ang 
# Date : 18/09/2019
# Proj : fast MNN pipeline

#' fast MNN pipeline
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
  }

  b_seurat <- FindVariableGenes(object = b_seurat, do.plot = TRUE, mean.function = ExpMean, dispersion.function = LogVMR, 
                                x.low.cutoff = b_x_low_cutoff, x.high.cutoff = b_x_high_cutoff, y.cutoff = b_y_cutoff)
  print(paste("Num var genes: ", length(x = b_seurat@var.genes)))
  if(length(x = b_seurat@var.genes) < numVG) {
    print("The number of highly variable genes in data is small, can change the threshold to gain more HVG genes")
  }

  # Get the list of highly variable gens in 2 batches
  hvg_all <- rownames(x = head(x = b_seurat@hvg.info, n = 5000))
  length(hvg_all)

  meta_data <- b_seurat@meta.data
  data_filtered <- b_seurat@scale.data[hvg_all,]

  #split data into batches
  data_list <- list()
  meta_list <- list()
  list_batch <- unique(meta_data[,batch_label])
  for (i in 1:length(list_batch)) {
    selected <- row.names(meta_data)[which (meta_data[, batch_label] == list_batch[i])]
    data_list[[i]] <- data_filtered[, selected]
    meta_list[[i]] <- meta_data[selected, ]
  }   

  return(list( data_list, meta_list ))
}

call_mnn <- function(data_list, meta_list, batch_label, celltype_label, npcs = 20, 
                         plotout_dir = "", saveout_dir = "", 
                         outfilename_prefix = "", 
                         visualize = T, save_obj = T)
{

  #plotting
  k_seed = 10
  tsne_perplex = 30
  tsneplot_filename = "_fastmnn_tsne"
  mnn_out_filename = "_fastmnn_out"
  metadata_out_filename = "_fastmnn_metadata_out"
  pca_filename = "_fastmnn_pca"

  b1_data <- data_list[[1]]
  b2_data <- data_list[[2]]
  b3_data <- data_list[[3]]

  b1_meta <- meta_list[[1]]
  b2_meta <- meta_list[[2]]
  b3_meta <- meta_list[[3]]

  metadata <- do.call(rbind, meta_list)

  ##########################################################
  #run

  t1 = Sys.time()

  cosd1 <- cosineNorm(as.matrix(b1_data))
  cosd2 <- cosineNorm(as.matrix(b2_data))
  cosd3 <- cosineNorm(as.matrix(b3_data))

  pcs_total <- multiBatchPCA(cosd1, cosd2, cosd3)   #defaut npca = 50
  out_mnn_total <- fastMNN(pcs_total[[1]], pcs_total[[2]], pcs_total[[3]], pc.input = TRUE)

  t2 = Sys.time()
  print(t2-t1)

  ##########################################################
  #save

  saveRDS(out_mnn_total, file=paste0(saveout_dir, outfilename_prefix, mnn_out_filename, ".RDS"))

  corrected_pca_mnn <- as.data.frame(out_mnn_total$corrected)
  rownames(corrected_pca_mnn) <- colnames(do.call(cbind, data_list))
  cells_use <- rownames(corrected_pca_mnn)
  metadata <- metadata[cells_use,]

  corrected_pca_mnn[rownames(metadata), batch_label] <- metadata[ , batch_label]
  #corrected_pca_mnn[rownames(metadata), 'batch'] <- metadata[ , "batch"]
  corrected_pca_mnn[rownames(metadata), celltype_label] <- metadata[ , celltype_label]
  write.table(corrected_pca_mnn, file=paste0(saveout_dir, outfilename_prefix, pca_filename, ".txt"), quote=F, sep='\t', row.names = T, col.names = NA)

  saveRDS(metadata, file = paste0(saveout_dir, outfilename_prefix, metadata_out_filename, ".RDS"))
  #write.table(metadata, file=paste0(saveout_dir, outfilename_prefix, metadata_out_filename, ".txt"), quote=F, sep="\t", row.name=TRUE, col.name=NA)

  if (visualize) {

    set.seed(10)
    out_tsne <- Rtsne(out_mnn_total$corrected, check_duplicates=FALSE, perplexity = 30)
    tsne_df<- as.data.frame(out_tsne$Y)

    rownames(tsne_df) <- colnames(do.call(cbind, data_list))
    dim(tsne_df)
    colnames(tsne_df) <- c('tSNE_1', 'tSNE_2')
    tsne_df$batchlb <- metadata[ , batch_label]
    #tsne_df$batch <- metadata[ , "batch"]
    tsne_df$CellType <- metadata[ , celltype_label]

    p01 <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2, colour = batchlb)) + geom_point(alpha = 0.6) + theme_bw()
    p01 <- p01 + labs(x='tSNE_1',y='tSNE_2',title='fast MNN')
    p01 <- p01 + theme(legend.title = element_text(size=17), 
                       legend.key.size = unit(1.1, "cm"),
                       legend.key.width = unit(0.5,"cm"), 
                       legend.text = element_text(size=14), 
                       plot.title = element_text(color="black", size=20, hjust = 0.5))

    p02 <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2, colour = CellType)) + geom_point(alpha = 0.6) + theme_bw()
    p02 <- p02 + labs(x='tSNE_1',y='tSNE_2',title='fast MNN')
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







