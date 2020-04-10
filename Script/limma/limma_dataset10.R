  library(limma)
  library(Seurat)
  packageVersion('Seurat')
  
  
  # Set working directory
  rm(list=ls())
  this.dir <- '/home/hoa/hoatran/demo_normalization/demoLimma/'
  setwd(this.dir)
  
  load('limma_utils.R') 

  # create base directory to save the results 
  base_name <- 'results_dataset10_hemato/'  # change to your dataset name
  dir.create(base_name, showWarnings = FALSE)

  data_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/dataset/dataset10_hematoMNN_Hoa/final/'

  sample_b1 <- read.table(paste0(data_dir,'b1_celltype.txt'), sep = "\t", row.names = 1, head=T, check.names = FALSE) ; dim(sample_b1)
  sample_b1$batch <- 'batch1'
  sample_b2 <- read.table(paste0(data_dir,'b2_celltype.txt'), sep = "\t", row.names = 1, head=T, check.names = FALSE) ; dim(sample_b2)
  sample_b2$batch <- 'batch2'

  data_b1 <- read.table(paste0(data_dir,"b1_exprs.txt"), sep = "\t", row.names = 1, head=T, check.names = FALSE) ; dim(data_b1)
  data_b2 <- read.table(paste0(data_dir,"b2_exprs.txt"), sep = "\t", row.names = 1, head=T, check.names = FALSE) ; dim(data_b2)

  # Verification
  sum(rownames(data_b1)==rownames(data_b2))
  sum(colnames(sample_b1)==colnames(sample_b2))
  myData <- cbind(data_b1, data_b2)
  mySample <- rbind(sample_b1, sample_b2)
  mySample$celltype <- mySample$CellType  
  mySample$batchlb  <- mySample$batch 

  # Filter data, return a Seurat object
  # already filtered cells and genes
  min_cells = 0  
  min_genes = 0 
  myFilteredData <- limma_filter_data(myData, mySample, min_cells, min_genes,regressUMI=FALSE, saveFilteredResult=TRUE)
  mySample <- mySample[colnames(myFilteredData),]
  dim(mySample)
  dim(myFilteredData)
  
  limma_srt <- limma_batch_effect_removal(myFilteredData, mySample, base_name, saveRDA=TRUE, savetxt=FALSE)
  

  limma_visualization(limma_srt, base_name, plotUMAP=FALSE, plotTSNE=TRUE, saveRDA=TRUE)
  
  