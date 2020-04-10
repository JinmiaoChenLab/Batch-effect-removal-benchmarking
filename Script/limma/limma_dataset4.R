  # Limma batch effect removal 
  # Input: log normalize matrix 
  # Documentation: https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/removeBatchEffect
  # Hoa Tran
  
  
  library(limma)
  library(Seurat)
  packageVersion('Seurat')
  
  # Set working directory
  rm(list=ls())
  this.dir <- '/home/hoa/hoatran/demo_normalization/demoLimma/'
  setwd(this.dir)
  
  load('limma_utils.R')  
  
  
  # create base directory to save the results 
  base_name <- 'results_dataset4_pancreatic/'  # change to your dataset name
  dir.create(base_name, showWarnings = FALSE)
  
  data_dir <- '/home/hoa/hoatran/demo_normalization/dataset/dataset4_human_pancreatic/filter_data_Hoa/'
  # source_dir <- '/home/hoa/hoatran/demo_normalization/source_code/evaluation/'
  # load(paste0(source_dir,'evaluation_utils.R'))  # change directory to the file
  # load(paste0(source_dir,'preprocess_utils.R'))  # change directory to the file
  
  TPM_file <- 'myData_pancreatic_5batches.txt'  # replace by link to dataset
  sample_file <- 'mySample_pancreatic_5batches.txt' # replace by link to dataset
  myData <- read.table(paste0(data_dir,TPM_file),sep="\t",header=T,row.names=1,check.names=F)
  mySample <- read.table(paste0(data_dir,sample_file),sep="\t",header=T,row.names=1,check.names=F)
 
  unique(mySample$batch)
  
  # Filter data, return a Seurat object
  min_cells = 10
  min_genes = 300 
  output_ls <- limma_filter_data(myData, mySample, min_cells, min_genes)
 
  
  myFilteredData <-  output_ls$pbmc@data[output_ls$hvg_genes,]
  dim(myFilteredData)
  mySample <- mySample[colnames(myFilteredData),]
  dim(mySample)
  
  limma_srt <- limma_batch_effect_removal(myFilteredData, mySample, base_name, saveRDA=TRUE, savetxt=TRUE)
  
  limma_visualization(limma_srt, base_name, plotUMAP=TRUE, plotTSNE=TRUE, saveRDA=TRUE)
  
  pbmc <- limma_srt
  