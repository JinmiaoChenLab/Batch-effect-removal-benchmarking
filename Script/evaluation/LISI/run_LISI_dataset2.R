library(ggplot2)
library(lisi)
library(ggpubr)

# Script for visualization 
rm(list=ls())

utils_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/source_code_github/evaluation_utils/'
source(paste0(utils_dir,'lisi/lisi_utils.R'))

this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset2_cellatlas/'
setwd(this_dir)
eval_metric <- 'LISI_v2/' 
dataset_use <- 'dataset2'
plx = 40
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
dir.create(paste0(this_dir, eval_metric), showWarnings = F)

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
# View(head(resnet_df))
# # Important: With resnet, need to check the cell name 
# resnet_df$cell <- gsub('-[0-9]$','',resnet_df$cell)
# rownames(resnet_df) <- resnet_df$cell


fn <- paste0(dataset_use,'_raw_pca.csv')
meta_ls <- get_celltype_common(data_dir, fn)
length(meta_ls$cells_common)

get_cells_integration_iLISI_v2(dataset_use, meta_ls, this_dir, plx, eval_metric)



# #######################################
######### cLISI
# #######################################

get_celltype_mixing_cLISI(dataset_use, this_dir, plx, eval_metric)

############################
## Summary output
############################


summary_LISI(this_dir, plottitle=paste0('LISI - ',dataset_use), plx, eval_metric)


# fn <- 'dataset2_raw_pca.csv'
# meta_ls <- get_celltype_common(data_dir, fn)
# length(meta_ls$ct_common)


# summary_LISI <- function(meta_ls, this_dir, plottitle='LISI - dataset', plx=40, eval_metric='LISI/',ht=400, wd=400){
#   iLISI_df <- read.csv(paste0(this_dir, eval_metric,"result/",plx,"/iLISI_summary.csv"), head=T, check.names = F)
#   colns <- c('methods_use', 'iLISI_median', 'iLISI_median_norm')
#   
#   median_iLISI <- normalize_values(iLISI_df, colns, min_val=min(iLISI_df), max_val=max(iLISI_df))
#   mini <- min(median_iLISI$iLISI_median_norm)
#   maxi <- max(median_iLISI$iLISI_median_norm)
#   median_iLISI$iLISI_median_norm2 <- (median_iLISI$iLISI_median_norm - mini) / (maxi - mini)
#   
#   
#   cLISI_df <- read.csv(paste0(this_dir,eval_metric,"result/",plx,"/cLISI_summary.csv"), head=T, check.names = F)
#   colns <- c('methods_use', 'cLISI_median', 'cLISI_median_norm')
#   median_cLISI <- normalize_values(cLISI_df, colns, min_val=min(cLISI_df), max_val=max(cLISI_df))
#   median_cLISI$cLISI_median_norm_sub <- 1 - median_cLISI$cLISI_median_norm
#   minc <- min(median_cLISI$cLISI_median_norm_sub)
#   maxc <- max(median_cLISI$cLISI_median_norm_sub)
#   median_cLISI$cLISI_median_norm_sub2 <- (median_cLISI$cLISI_median_norm_sub - minc) / (maxc - minc)
#   
#   
#   final_df = merge(median_iLISI, median_cLISI, by="methods_use")
#   final_df$sum_normXY <- final_df$iLISI_median_norm2 + final_df$cLISI_median_norm_sub2
#   final_df$fscore <- (2 * final_df$iLISI_median_norm2 * final_df$cLISI_median_norm_sub2)/
#                      (final_df$iLISI_median_norm2 + final_df$cLISI_median_norm_sub2)
#   
#   final_df <- final_df[order(final_df$fscore, decreasing = T),]
#   # final_df$cLISI_median_norm_sub = 1 - final_df$cLISI_median_norm
#   write.csv(final_df,paste0(this_dir, eval_metric, "result/", plx, '/', 'summary_median_', plx, '.csv'), 
#               quote=F, row.names=F)
#   
# 
#   # plot final LISI
#   plot_final_LISI(final_df, plottitle, this_dir, plx, ht, wd, eval_metric, 
#                   xstring = 'cLISI_median_norm_sub', ystring = 'iLISI_median_norm', plottype = 'methods_use')
#   
# }

