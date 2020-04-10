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
    write.table(myRunTime, file = paste0(base_name,labeluc,"_runtime.txt"), row.names = TRUE, col.names = TRUE, sep="\t")
  }  

  # can use 5000 genes here 
  # LIMMA : removeBatchEffect(x, batch=NULL, batch2=NULL, covariates=NULL, design=matrix(1,ncol(x),1), ...)
  # x: numeric matrix, or any data object that can be processed by getEAWP containing log-expression values
  # for a series of samples. Rows correspond to probes and columns to samples.
  # batch: factor or vector indicating batches.
  # myData: contain log-expression values
  # Hoa Tran
  
  limma_batch_effect_removal <- function(myData, mySample, base_name, saveSeuratObjRDS=FALSE, save_mtx=TRUE){
    
    library(ezTools)
    mySample <- mySample[colnames(myData),]
    t1 = Sys.time()
    lm_df <- removeBatchEffect(as.matrix(myData), factor(mySample$batch))
    t2 = Sys.time()
    print(t2-t1)
    
    
    ######################
    ####  Execution time 
    ######################
    # Export execution time 
    labeluc <- "limma_time_execution"
    runtime_export_func(t1, t2, labeluc, base_name)
    
    # Save output
    if(save_mtx){
      fn <- "limma_output_matrix.rds"
      saveRDS(lm_df, file = paste0(base_name,fn))
      # write.table(lm_df, file=paste0(base_name,"limma_output.txt"), sep="\t", row.name=TRUE, col.name=TRUE) 
      ezTools::fast_save_table(lm_df, base_name, fn, row.names = T, col.names = T, quote = F, sep = "\t")
      return(lm_df)
    }
    
    
    if(saveSeuratObjRDS){
      limma_srt <- CreateSeuratObject(raw.data = lm_df, meta.data = mySample, project = "Limma_normalized")
      saveRDS(limma_srt, file = paste0(base_name,"limma_normalized_seurat_obj.rds"))
      return(limma_srt)
    }
    
    
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
  
  limma_visualization <- function(pbmc, base_name='',plotUMAP=TRUE, plotTSNE=TRUE, saveRDS=TRUE){
    # pbmc@data <- pbmc@raw.data
    pbmc@scale.data <- pbmc@data
   
    npcs <- 20
    pbmc <- RunPCA(object = pbmc, pcs.compute = npcs, pc.genes = rownames(pbmc@data), do.print = TRUE, pcs.print = 1:5, 
                   genes.print = 5) 
    if(plotUMAP){
      pbmc <- RunUMAP(pbmc, reduction.use = "pca", dims.use = 1:20)  
      png(paste(base_name,"umap_limma.png",sep=""),width = 2*800, height = 2*500, res = 2*72)
      p11 <- DimPlot(object = pbmc, reduction.use = 'umap', group.by ="batchlb")
      # p12 <- DimPlot(object = pbmc, reduction.use = 'umap', group.by ="celltype")   # or orig.ident
      print(p11)
      dev.off()
    } 
    if(plotTSNE){
      pbmc <- RunTSNE(pbmc, reduction.use = "pca", dims.use = 1:20, do.fast = T, perplexity=30, check_duplicates=F)
      png(paste(base_name,"tsne_limma.png",sep=""),width = 2*1100, height = 2*800, res = 2*72)
      p1 <- TSNEPlot(pbmc, do.return = T, pt.size = 0.5, group.by = "batchlb")
      p2 <- TSNEPlot(pbmc, do.return = T, pt.size = 0.5, group.by = "celltype")
      print(p1)
      dev.off()
    }
    if(saveRDS){
      saveRDS(pbmc, file = paste0(base_name,"limma_normalized_seurat_obj_umap_tsne.rda"))
    }
    
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
  
  