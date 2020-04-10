
library(ggplot2)
library(lisi)
library(ggpubr)

# Script for visualization 
rm(list=ls())
this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset1_uc3/'

setwd(this_dir)
eval_metric <- 'LISI_v2/' 

utils_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/source_code_github/evaluation_utils/'
source(paste0(utils_dir,'lisi/lisi_utils.R'))

dataset_use <- 'dataset1'
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


fn <- 'dataset1_raw_pca.csv'
meta_ls <- get_celltype_common(data_dir, fn)
length(meta_ls$cells_common)
meta_ls$ct_common

get_cells_integration_iLISI_v2(dataset_use, meta_ls, this_dir, plx, eval_metric)


# #######################################
######### cLISI
# #######################################

get_celltype_mixing_cLISI(dataset_use, this_dir, plx, eval_metric)




summary_LISI(this_dir, plottitle=paste0('LISI - ',dataset_use), plx, eval_metric)







# Test different perplexities
# minb = 0
# maxb = length(unique(myPCA$batch))
# minct = 1
# maxct = length(unique(myPCA$cell_type))
# 
# median_batch_norm <- c()
# median_ct_norm <- c()
# plxls <- c()
# 
# for(plx in seq(10, 150, 10)) {
#   plxls <- c(plxls, plx)
#   print('Run LISI with: ')
#   print(plx)
#   lisi_res <- lisi::compute_lisi(lisi_embeddings, lisi_meta_data, lisi_label,perplexity = plx)
#   lisi_res$cell <- rownames(lisi_embeddings)
#   
#   lisi_batch <- subset(lisi_res,select=c('batch','cell'))
#   lisi_celltype <- subset(lisi_res,select=c('cell_type','cell'))
#   
#   median_batch_norm <- c(median_batch_norm, (median(lisi_batch[,'batch']) - minb) / (maxb - minb))
#   median_ct_norm <- c(median_ct_norm, 1-((median(lisi_celltype[,'cell_type']) - minct) / (maxct - minct)))
# }
# 
# LISI_df <- data.frame('perplexity_use'=plxls, 'method'=rep(methods_use, length(plxls)), 
#                       'iLISI_median_batch'=median_batch_norm,
#                       'cLISI_median_celltype'=median_ct_norm)
# write.table(LISI_df, paste0(save_dir,eval_metric,methods_use,'/','lisi_median_perplexities.txt'), 
#             quote=F, sep='\t', row.names=T, col.names=NA)
# return(LISI_df)

# finaldf = rbind(limma_df, harmony_df, seurat2_df, scgen_df, resnet_df, combat_df,
#                 scanorama_df, fastMNN_df, seurat3_df)
# 
# dim(finaldf)
# write.table(finaldf, paste0(this_dir,eval_metric,'/','lisi_median_perplexities_summary.txt'), 
#             quote=F, sep='\t', row.names=F, col.names=T)
# 
# dftmp <- finaldf[which(!finaldf$method %in% 'seurat2_multicca'),]
# dim(dftmp)
# pcLISI <- lisi_perplexity_plot(dftmp, xstring='method', 
#                                ystring='cLISI_median_celltype',
#                                plottype = 'perplexity_use',
#                                toptitle='LISI perplexity cell type', title='1 _cLISI cell type')
# 
# pcLISI
# 
# piLISI <- lisi_perplexity_plot(dftmp, xstring='method', 
#                                ystring='iLISI_median_batch',
#                                plottype = 'perplexity_use',
#                                toptitle='LISI perplexity - batch', title='iLISI cell type')
# 
# piLISI
# 
# 
# 
# library(gridExtra)
# lay <- rbind(c(1,2))
# 
# base_name = 'eval'
# plots <- list(piLISI, pcLISI)
# png(paste(this_dir,"LISI_perplexity.png",sep=""),height = 2*350, width=2*850,res = 2*72)
# # pdf(paste(data_dir,"dataset7_iLISI_cLISI_kBET.pdf",sep=""),height=plotheight, width=plotwidth)
# print(grid.arrange(grobs = plots, layout_matrix = lay,top=" ",
#                    bottom=" ",right=" "))
# dev.off()
