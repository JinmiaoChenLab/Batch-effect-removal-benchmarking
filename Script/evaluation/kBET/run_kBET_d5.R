# Run all kBET

library(ggplot2)
library(kBET)
library(ggpubr)
library(ezTools)


##################
## dataset 5
##################

rm(list=ls())
this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset5_human_pbmc/'
setwd(this_dir)
data_dir = '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/dataset5/'
eval_metric <- 'kBET_Hoa_v2/'
# kn = 30
dataset_use <- 'dataset5'

utils_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/source_code_github/evaluation_utils/'
source(paste0(utils_dir,'evaluation_utils.R'))
source(paste0(utils_dir,'kbet/kbet_utils.R'))


method_dir <- c('raw','seurat2','seurat3',
                'harmony','fastMNN','classicMNN',
                'combat','limma','scgen',
                'scanorama','resnet','zinbwave',
                'scmerge','liger')


methods_use <- c('Raw','Seurat_2','Seurat_3',
                 'Harmony','fastMNN','MNN_Correct',
                 'ComBat','limma','scGen',
                 'Scanorama','MMD-ResNet','ZINB-WaVE',
                 'scMerge','LIGER')


# extract cells use
dir.create(paste0(this_dir, eval_metric), showWarnings = F)
fn <- paste0(dataset_use, '_raw_pca.csv')
meta_ls <- get_celltype_common_kbet_V2(data_dir, fn)


#### iteration 
size = length(meta_ls$cells_extract)  # 13908
pct = c(seq(from = 5, to = 25, by = 5))
kns = floor(pct*size/100)


res_dir <- paste0(this_dir, eval_metric,"result/")


###################
## kbet computation
###################
library(parallel)
mclapply(1:length(kns), function(i){
  kn = kns[i]
  print(paste0('____k0 value use is: ', kn))
  summary_KBET_v2(meta_ls, data_dir, methods_use, method_dir, dataset_use, eval_metric, this_dir, kn)
}, mc.cores = 6)


###################
## kbet generate results
###################

generate_median_kbet(pct, kns, res_dir, dataset_use, this_dir, kbet_output=T)


# pct = c(2:10, seq(from = 15, to = 50, by = 5))
# kns = floor(pct * size/100)
# uc <- paste0(pct,"%")
# for(i in 1:length(kns)){
#   kn = kns[i]
#   pkbet <- get_kbet_output(dataset_use, this_dir, kn, eval_metric, save_results=T)
#   
#   dir_result <- paste0(res_dir, kn, "/Acceptance rate_median_val_kBET.csv")
#   
#   if(file.exists(dir_result)){
#     res_i = read.csv(dir_result, header=T)
#     print(colnames(res_i))
#     rej_i <- subset(res_i, select = c("methods","median_observed"))
#     acc_i <- subset(res_i, select = c("methods","median_acceptance_rate"))
#     print(dim(rej_i))
#     if(i==1){
#       rej = rej_i
#       acc = acc_i
#     }else{
#       rej = merge(rej, rej_i, by="methods")  
#       acc = merge(acc, acc_i, by="methods") 
#     }
#   }  
# }
# 
# colnames(rej) <- c("methods", uc)
# colnames(acc) <- c("methods", uc)
# med_rej <- c()
# med_acc <- c()
# for(r in 1:nrow(rej)){
#   med_rej <- c(med_rej, median(as.numeric(rej[r,uc])))
#   med_acc <- c(med_acc, median(as.numeric(acc[r,uc])))
# }
# rej$median_val <- med_rej
# acc$median_val <- med_acc
# 
# rej = rej[order(rej$median_val, decreasing = F),]
# acc = acc[order(acc$median_val, decreasing = T),]
# 
# write.csv(rej, paste0(res_dir,dataset_use,"_kbet_rejection.csv"))
# write.csv(acc,paste0(res_dir,dataset_use,"_kbet_acception.csv"))
