library(ggplot2)
library(lisi)
library(ggpubr)

# Script for visualization 
rm(list=ls())

utils_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/source_code_github/evaluation_utils/'
source(paste0(utils_dir,'lisi/lisi_utils.R'))


this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset5_human_pbmc/'
setwd(this_dir)
eval_metric <- 'LISI_v2/' 
dataset_use <- 'dataset5'
plx = 40 
dir.create(paste0(this_dir, eval_metric), showWarnings = F)
# Get output of LISI
data_dir = paste0('/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/',dataset_use,'/')
methods_use <- c('Raw','Seurat_2','Seurat_3',
                 'Harmony','fastMNN','MNN_Correct',
                 'ComBat','limma','scGen',
                 'Scanorama','MMD-ResNet','ZINB-WaVE',
                 'scMerge','LIGER')

method_dir <- c('raw','seurat2','seurat3',
                'harmony','fastMNN','classicMNN',
                'combat','limma','scgen',
                'scanorama','resnet','zinbwave',
                'scmerge','liger')

fn_ls <- c()
for(i in rep(1:length(method_dir),1)){
  fn_ls <- c(fn_ls, paste0(dataset_use,'_', method_dir[i],'_','pca.csv'))
}
print(fn_ls)



length(methods_use) == length(fn_ls)
for(i in rep(1:length(methods_use), 1)){
  print(methods_use[i])
  run_LISI_final(fn_ls[i], data_dir, this_dir, eval_metric, methods_use[i], plx)
}


#####################################
# iLISI batch
# Combine all results together
# 13 methods (except BBKNN) + raw data
#####################################



# rownames(resnet_df)[1:5]
# Important: With resnet, need to check the cell name 
# For this dataset, we do not do this
# resnet_df$cell <- gsub('-[0-9]$','',resnet_df$cell)
# rownames(resnet_df) <- resnet_df$cell
# print(rownames(resnet_df)[1:3])
# print(resnet_df$cell[1:3])

fn <- 'dataset5_raw_pca.csv'
data_dir = '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/dataset5/'
meta_ls <- get_celltype_common(data_dir, fn)
length(meta_ls$ct_common)


get_cells_integration_iLISI_v2(dataset_use, meta_ls, this_dir, plx, eval_metric)



# #######################################
######### cLISI
# #######################################


# Important: With resnet, need to check the cell name 
# resnet_df$cell <- gsub('-[0-9]$','',resnet_df$cell)
# rownames(resnet_df) <- resnet_df$cell
# print(rownames(resnet_df)[1:3])
# print(resnet_df$cell[1:3])


get_celltype_mixing_cLISI(dataset_use, this_dir, plx, eval_metric)





############################
## Summary output
############################
summary_LISI(this_dir, plottitle=paste0('LISI - ',dataset_use), plx, eval_metric)
