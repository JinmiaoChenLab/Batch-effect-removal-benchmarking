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
    myData=myData[,cells_use]
  }
  if(is_filter_genes){
    # Gene filtering
    num_cells = rowSums(myData > 0)            
    genes_use = names(num_cells[which(num_cells>min_cells)])
    myData = myData[genes_use,]
  }
  
  # Normalize each library to the median of the transcript counts across all cells 
  # Then, log transform expression values   
  print("Median normalizing counts and log-transforming")
  col_sums = apply(myData,2, sum)
  med_trans = median(col_sums)
  norm_counts = med_trans* scale(myData, center=FALSE, scale=col_sums)
  rm(myData)
  myFilteredData = as.data.frame(log(norm_counts + 1))
  
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
                        scale.factor = 100)
  
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
## COMBAT function
#########################################

# Opt1: parametric option: mod=NULL, par.prior=TRUE, prior.plots=FALSE, mean_only=FALSE  the best option
# Opt2: non-parametric adjustment, mean-only version mod=NULL, par.prior=FALSE,prior.plots=FALSE, mean.only=TRUE
# Pls load evaluation functions first to have time_execution function

# myData: data frame or matrix, genes x cells 
# batch: a vector of bach label, ex: batch <- pbmc@meta.data[cells_use,'batch']
# In case big data, can save as rda or rds format, save_txt=FALSE
combat_batch_effect_removal <- function(myData, batch=NULL, mod=NULL, par_prior=TRUE, prior_plots=FALSE, mean_only=FALSE, 
                                        base_name='parametric/', save_txt=FALSE, saveRDS=TRUE){
  library("sva")
  # Create a folder to save the results
  dir.create(base_name, showWarnings = FALSE)
  
  
  # Run COMBAT
  # Pls run 2 options: parametric and non parametric and get the best output
  t1 = Sys.time()
  combat_output = ComBat(dat=as.matrix(myData), batch=batch, 
                         mod=mod, par.prior=par_prior, mean.only=mean_only, prior.plots=prior_plots)
  t2 = Sys.time()
  
  if(save_txt){
    write.table(combat_output, file = paste0(base_name,"combat_normalized_matrix.txt"), row.names = T, col.names = T, sep="\t")  
    print("Data saved as txt format")
    
    # library(ezTools)
    # fn = 'combat_normalized_matrix.txt'
    # ezTools::fast_save_table(combat_output, base_name, fn, 
    #                          row.names = T, col.names = T, quote = F, sep = "\t")
  }
  if(saveRDS){
    saveRDS(combat_output, file = paste0(base_name,"combat_normalized.rds"))
    print("Data saved as rda format")
  }
  
  
  ######################
  ####  Execution time 
  ######################
  # Export execution time, pls load this function from evaluation_utils.R first 
  labeluc <- "combat_time_execution"
  runtime_export_func(t1, t2, labeluc, base_name)
  return(combat_output)
}
