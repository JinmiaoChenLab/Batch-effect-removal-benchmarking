library(ggplot2)
library(lisi)
library(ggpubr)
library(stringr)

rm(list=ls())



this_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/manuscript_results/dataset6_cell_line/'
setwd(this_dir)
# eval_metric <- 'LISI/' 
eval_metric <- 'LISI_v2/' 
dataset_use <- 'dataset6'
plx = 40 

utils_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/source_code_github/evaluation_utils/'
source(paste0(utils_dir,'lisi/lisi_utils.R'))

# Get output of LISI
data_dir = '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/dataset6/'
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
  fn_ls <- c(fn_ls, paste0(dataset_use,'_', method_dir[i],'_pca.csv'))
}
print(fn_ls)
dir.create(paste0(this_dir, eval_metric), showWarnings = F)

library(parallel)
mclapply(1:length(methods_use), function(i){
# for(i in rep(1:length(methods_use), 1)){
  print(methods_use[i])
  run_LISI_final_d6(fn_ls[i], data_dir, this_dir, eval_metric, methods_use[i], c(1,3), blabel='13', plx)
  run_LISI_final_d6(fn_ls[i], data_dir, this_dir, eval_metric, methods_use[i], c(2,3), blabel='23', plx)
}, mc.cores = 12)




mclapply(1:length(methods_use), function(i){
# for(i in rep(1:length(methods_use), 1)){
  print(methods_use[i])
  run_LISI_final_celltype(fn_ls[i], data_dir, this_dir, eval_metric, methods_use[i], plx)
}, mc.cores = 12)


#####################################
# iLISI batch
# Combine all results together
# 13 methods (except BBKNN) + raw data
#####################################

get_cells_integration_iLISI_d6(dataset_use, '13', this_dir, plx, eval_metric)
get_cells_integration_iLISI_d6(dataset_use, '23', this_dir, plx, eval_metric)



# #######################################
######### cLISI
# #######################################

# #######################################
######### cLISI
# #######################################



get_celltype_mixing_cLISI(dataset_use, this_dir, plx, eval_metric)



############################
## Summary output
############################


summary_LISI_d6(dataset_use, this_dir, plx, eval_metric)


