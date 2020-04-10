library(ggplot2)
library(lisi)
library(ggpubr)
library(ezTools)


rm(list=ls())

this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset9_Human_cell_atlas/'
setwd(this_dir)
eval_metric <- 'LISI_v2/' 
dataset_use <- 'dataset9'

utils_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/source_code_github/evaluation_utils/'
source(paste0(utils_dir,'lisi/lisi_utils.R'))

data_dir = paste0('/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/',dataset_use,'/')

plx = 40 
dir.create(paste0(this_dir, eval_metric), showWarnings = F)

############################
## LISI computation
############################

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


library(parallel)
mclapply(1:length(methods_use), function(i){
# for(i in rep(1:length(methods_use), 1)){
  print(methods_use[i])
  run_LISI_final_dataset9(fn_ls[i], data_dir, this_dir, eval_metric, methods_use[i], plx)
}, mc.cores = 12)




#####################################
# iLISI batch
# Combine all results together
# 13 methods (except BBKNN) + raw data
#####################################
# this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset8_Mouse_brain/'
# setwd(this_dir)


seurat2_df <- read.table(paste0(this_dir,eval_metric,"Seurat_2/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
seurat3_df <- read.table(paste0(eval_metric,"Seurat_3/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
harmony_df <- read.table(paste0(eval_metric,"Harmony/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
fastMNN_df <- read.table(paste0(eval_metric,"fastMNN/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
resnet_df <- read.table(paste0(eval_metric,"MMD-ResNet/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
scanorama_df <- read.table(paste0(eval_metric,"Scanorama/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
scGen_df <- read.table(paste0(eval_metric,"scGen/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
raw_df <- read.table(paste0(eval_metric,"Raw/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
correctMNN_df <- read.table(paste0(eval_metric,"MNN_Correct/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
combat_df <- read.table(paste0(eval_metric,"ComBat/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
liger_df <- read.table(paste0(eval_metric,"LIGER/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
limma_df <- read.table(paste0(eval_metric,"limma/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
scMerge_df <- read.table(paste0(eval_metric,"scMerge/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)
zinbwave_df <- read.table(paste0(eval_metric,"ZINB-WaVE/lisi_batch","_",plx,".txt"), head=T, row.names = 1, check.names = FALSE)


dir.create(paste0(this_dir, eval_metric,"result/"), showWarnings = F)
dir.create(paste0(this_dir, eval_metric,"result/",plx,"/"), showWarnings = F)
methods_use <- c('Raw','Seurat_2','Seurat_3','Harmony','fastMNN','MNN_Correct','ComBat',
                 'limma','scGen','Scanorama','MMD-ResNet','ZINB-WaVE','scMerge','LIGER')
piLISI <- LISI_boxplot_fun(data = list(raw_df, seurat2_df, seurat3_df,
                                       harmony_df, fastMNN_df, correctMNN_df,
                                       combat_df, limma_df, scGen_df,
                                       scanorama_df, resnet_df,
                                       zinbwave_df, scMerge_df, liger_df),
                           vect_names_method = methods_use,
                           title='iLISI', toptitle=paste0('iLISI Batch Mixing plx ',plx),
                           base_name=paste0(this_dir, eval_metric, "result/",plx,"/"), save_results=T)
piLISI


save(piLISI, file = paste0(this_dir, eval_metric,"result/",plx, "/piLISI_dataset8_plots_plx_",plx,".rda"))




summary_LISI_d9(this_dir, plottitle=paste0('LISI - ',dataset_use), plx, eval_metric)



iLISI_df <- read.csv(paste0(this_dir, eval_metric,"result/",plx,"/iLISI_summary.csv"), head=T, check.names = F)
min_val=min(iLISI_df)
max_val=max(iLISI_df)
norm_iLISI_df <- (iLISI_df - min_val) / (max_val - min_val)
write.csv(norm_iLISI_df, paste0(this_dir, eval_metric, "result/", plx, '/', 'summary_norm_median_', plx, '.csv'),
          quote=F, row.names=F)



