rm(list=ls())

library(limma)
library(Seurat)  # Seurat 2 version
packageVersion('Seurat')
library(ezTools)
  
source('limma_utils.R')

# Set working directory
  
this.dir <- '/home/chenjm/Hoa_batch_normalization/demo_limma/'
setwd(this.dir)
  
# create base directory to save the results 
base_name <- 'results_dataset9_HCA/'  # change to your dataset name
dir.create(base_name, showWarnings = FALSE)
  
data_dir <- '/home/xm/Projects/batch_norm/datasets/dataset9_Human_cell_atlas/filtered_genes_and_cells/'
TPM_file <- 'HCA_genes_cells_filtered_filtered_UMI.RDS'  # replace by link to dataset
sample_file <- 'HCA_genes_cells_filtered_filtered_cell_info.txt' # replace by link to dataset
  
mySample <- read.table(file = paste0(data_dir,sample_file),sep="\t",header=T,row.names=1,check.names = F)

mySample$batch <- ifelse(mySample$tissue=="Bone marrow","batch1",
                         ifelse(mySample$tissue=="Cord blood","batch2",NA))



summary(as.factor(mySample$batch))
  
# Process sample file, done
# ext <- c("donor","tissue")
# mySample <- mySample[ , (names(mySample) %in% ext)]
# mySample$batch <- ifelse(mySample$tissue=="Bone marrow",1,
#                          ifelse(mySample$tissue=="Cord blood",2,NA))
# 
# mySample$batchlb <- ifelse(mySample$tissue=="Bone marrow","Bone_Marrow",
#                          ifelse(mySample$tissue=="Cord blood","Cord_Blood",NA))
# 
# unique(mySample$batchlb)
# # mySample$celltype <- 'NA'
# write.table(mySample, file = paste0(base_name,"HCA_genes_cells_filtered_filtered_cell_info_correct.txt"), row.names = TRUE, col.names = TRUE, sep="\t")
# 
# percent_extract <- 0.2
# savefnRDS <- paste0('downsample_HCA_',percent_extract,".rds")
  
    
# hvg_genes_fn <- 'HCA_variable_genes.txt'
# myHVGgenes <- read.table(file = paste0(data_dir,hvg_genes_fn),sep="\t",header=T,row.names=1,check.names = F)
# head(myHVGgenes)
# dim(myHVGgenes)
# hvg_genes <- myHVGgenes[,1]
# length(hvg_genes)
# head(hvg_genes)

  
myData <- readRDS(paste0(data_dir,TPM_file))
dim(myData)
# min_cells = 100  # XM hv already filtered
# min_genes = 0
  
##########################################
##### Filtering 
##########################################
# limma use log transform values as input
# Get 5000 HVG genes and run limma
myFilteredData <- limma_filter_data(myData, mySample)
dim(myFilteredData)
dim(mySample)
rm(myData)  
  

  ##########################################
  ##### LIMMA batch effect removal
  ##########################################
  limma_df <- limma_batch_effect_removal(myFilteredData, mySample, base_name,  saveSeuratObjRDS=FALSE, save_mtx=TRUE, bigdata=T)
  
  
  # Normal workflow for small dataset using Seurat v2
  # # First get seurat object
  # limma_srt <- limma_batch_effect_removal(myFilteredData, mySample, base_name, saveSeuratObjRDS=TRUE, save_mtx=FALSE)
  # # Then visualize seurat output
  # limma_visualization(limma_srt, base_name, plotUMAP=FALSE, plotTSNE=TRUE, saveRDA=TRUE)
  
 