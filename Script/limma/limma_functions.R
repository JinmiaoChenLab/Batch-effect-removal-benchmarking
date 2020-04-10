###### break sparse matrix into trunk, convert into dense matrix and then concatenate

preprocess_big <- function(sparse_matrix=norm_data,by=100000){
  n_col = ncol(sparse_matrix)
  n = floor(n_col/by)
  
  if(n<1){
    return(as.matrix(sparse_matrix))
  }else{
    for(i in 1:n){
      mat <- as.matrix(sparse_matrix[,((i-1)*by+1):(i*by)])
      
      if(i<2){
        res <- mat
      }else{
        res <- cbind(res,mat)
      }
      rm(mat)
    }
    if(n_col > n*by){
      res <- cbind(res,as.matrix(sparse_matrix[,(n*by+1):n_col]))
    }
  }
  rm(sparse_matrix)
  return(res)
}

#########################################
# Filter data
#########################################

filter_data_mtx <- function(myData, base_name, is_filter_cells=FALSE,min_genes=300, 
                            is_filter_genes=FALSE, min_cells=10){
  
  if(is_filter_cells){
    # Cell filtering
    num_genes = colSums(myData > 0); 
    names(num_genes) = colnames(myData)
    cells_use = names(num_genes[which(num_genes>min_genes)])
    myFilteredData=myData[,cells_use]
  }
  if(is_filter_genes){
    # Gene filtering
    num_cells = rowSums(myData > 0)            
    genes_use = names(num_cells[which(num_cells>min_cells)])
    myFilteredData = myData[genes_use,]
  }
  
  return(myFilteredData)
  
}

filter_data_seurat <- function(myData, mySample, min_cells=10, min_genes=300, group_col='celltype', regressUMI=FALSE){
  
  
  
  pbmc <- CreateSeuratObject(raw.data = myData, project = "Harmony_benchmark", min.cells = min_cells, min.genes = min_genes)
  dim(pbmc@data)
  
  
  cells_use <- colnames(pbmc@raw.data)
  length(cells_use)
  dim(mySample)
  # mySample_backup <- mySample
  mySample <- mySample[cells_use,]
  pbmc@meta.data$batch <- mySample$batch  
  # pbmc@meta.data$batchlb <-  paste0('Batch_',pbmc@meta.data$batch)
  pbmc@meta.data$batchlb <- mySample$batchlb
  pbmc@meta.data$celltype <- mySample$celltype
  
  
  
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
  
  if(!regressUMI){  # by default
    pbmc <- ScaleData(object = pbmc)
  }else{
    pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))   # in case read count, need to normalize by sum of each column (CPM normalization)
  }
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  # if(length(x = pbmc@var.genes) < 400)
  # {
  #   print("The number of highly variable genes is small, can change the threshold to gain more HVG genes")
  # }  
  
  return(pbmc)
  
}

limma_filter_data <- function(myData, mySample, min_cells=10, min_genes=300, regressUMI=FALSE, saveFilteredResult=FALSE){
  
  cells_use <- colnames(myData)
  mySample <- mySample[cells_use,]
  pbmc <- CreateSeuratObject(raw.data = myData, meta.data = mySample, project = "Limma_benchmark", min.cells = min_cells, min.genes = min_genes)
  
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                        scale.factor = 1e4)
  
  if(!regressUMI){  # by default
    pbmc <- ScaleData(object = pbmc)
  }else{
    pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))   # in case read count, need to normalize by sum of each column (CPM normalization)
  }
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.1)
  # if(length(x = pbmc@var.genes)<400)
  # {
  #   print("The number of highly variable genes is small, can change the threshold to gain more HVG genes")
  # }  
  
  nb_hvg <- 5000  # select only 5000 hvg genes
  if(nb_hvg>nrow(pbmc@data)){
    nb_hvg <- nrow(pbmc@data)
  }
  hvg_genes <- rownames(x = head(x = pbmc@hvg.info, n = nb_hvg))
  print(length(hvg_genes))
  if(saveFilteredResult){
    saveRDS(list(pbmc,hvg_genes), "filtered_seurat_obj.rds")
  }
  output_df <- pbmc@data[hvg_genes,]   # limma use log transform values as input
  
  
  return(output_df)
  
}
#########################################
# time util function
#########################################

runtime_export_func <- function(t1, t2, labeluc, base_name, timeunits='secs_mins')
{
  label_usecase = labeluc  #label "usecase3_method_abc"
  time_secs = c()
  time_secs = c(time_secs, as.numeric(difftime(t2,t1,units='secs')))
  print(tail(time_secs,n=1))
  
  time_mins = c()
  time_mins = c(time_mins, as.numeric(difftime(t2,t1,units='mins')))
  print(tail(time_mins,n=1))
  
  myRunTime <- data.frame("use_case"=label_usecase, "exetime_secs"=time_secs, "exetime_mins"=time_mins)
  write.table(myRunTime, file = paste0(base_name,labeluc,"_runtime.txt"), row.names = FALSE, col.names = TRUE, sep="\t")
}  



#########################################
## limma function
#########################################
# can use 5000 genes here 
# LIMMA : removeBatchEffect(x, batch=NULL, batch2=NULL, covariates=NULL, design=matrix(1,ncol(x),1), ...)
# x: numeric matrix, or any data object that can be processed by getEAWP containing log-expression values
# for a series of samples. Rows correspond to probes and columns to samples.
# batch: factor or vector indicating batches.
# myData: contain log-expression values
# In case big data, can save as rda or rds format, save_txt=FALSE

limma_batch_effect_removal <- function(myData, mySample, base_name, save_txt=FALSE, saveRDS=TRUE){
  
  # Create a folder to save the results
  dir.create(base_name, showWarnings = FALSE)
  
  mySample <- mySample[colnames(myData),]
  t1 = Sys.time()
  limma_srt <- removeBatchEffect(as.matrix(myData), factor(mySample$batch))
  t2 = Sys.time()
  
  if(save_txt){
    write.table(limma_srt, file = paste0(base_name,"limma_normalized_matrix.txt"), row.names = T, col.names = T, sep="\t")  
    print("Data saved as txt format")
    
    # library(ezTools)
    # fn = 'combat_normalized_matrix.txt'
    # ezTools::fast_save_table(combat_output, base_name, fn, 
    #                          row.names = T, col.names = T, quote = F, sep = "\t")
  }
  if(saveRDS){
    saveRDS(limma_srt, file = paste0(base_name,"limma_normalized.rds"))
    print("Data saved as rda format")
  }
  
  
  ######################
  ####  Execution time 
  ######################
  # Export execution time, pls load this function from evaluation_utils.R first 
  labeluc <- "limma_time_execution"
  runtime_export_func(t1, t2, labeluc, base_name)
  return(limma_srt)
}

limma_visualization <- function(pbmc, base_name='',plotUMAP=TRUE, plotTSNE=TRUE, saveRDS=TRUE){
  pbmc@data <- pbmc@raw.data
  pbmc@scale.data <- pbmc@raw.data
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.1)
  length(pbmc@var.genes)
  npcs <- 20
  
  pbmc <- RunPCA(object = pbmc, pcs.compute = npcs, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
                 genes.print = 5) 
  if(plotUMAP){
    pbmc <- RunUMAP(pbmc, reduction.use = "pca", dims.use = 1:20)  
    png(paste(base_name,"umap_limma.png",sep=""),width = 2*800, height = 2*500, res = 2*72)
    p11 <- DimPlot(object = pbmc, reduction.use = 'umap', group.by ="batchlb")
    p12 <- DimPlot(object = pbmc, reduction.use = 'umap', group.by ="celltype")   # or orig.ident
    print(plot_grid(p11, p12))
    dev.off()
  } 
  if(plotTSNE){
    pbmc <- RunTSNE(pbmc, reduction.use = "pca", dims.use = 1:20, do.fast = T, perplexity=30)
    png(paste(base_name,"tsne_limma.png",sep=""),width = 2*1100, height = 2*800, res = 2*72)
    p1 <- TSNEPlot(pbmc, do.return = T, pt.size = 0.5, group.by = "batchlb")
    p2 <- TSNEPlot(pbmc, do.return = T, pt.size = 0.5, group.by = "celltype")
    print(plot_grid(p1, p2))
    dev.off()
  }
  if(saveRDS){
    saveRDS(pbmc, file = paste0(base_name,"limma_normalized_seurat_obj_umap_tsne.rda"))
  }
  
}
