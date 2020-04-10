library(ggplot2)
library(kBET)
library(ggpubr)
library(ezTools)


##################
## dataset 6
##################

rm(list=ls())
this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset6_cell_line/'
setwd(this_dir)
data_dir = '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/dataset6/'
eval_metric <- 'kBET_Hoa_v2/' 
# kn = 30
dataset_use <- 'dataset6'

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

dir.create(paste0(this_dir, eval_metric), showWarnings = F)
fn <- paste0(dataset_use, '_raw_pca.csv')
meta_ls_13 <- get_common_celltype_d6(c(1,3), data_dir, fn)
meta_ls_23 <- get_common_celltype_d6(c(2,3), data_dir, fn)

# pct = c(2:10)
# pct = c(5)
pct = c(seq(from = 10, to = 25, by = 5))
mean_sample_size = (length(meta_ls_13$cells_extract)+length(meta_ls_23$cells_extract))/2  # 3388
kns = floor(pct *  mean_sample_size / 100)


res_dir <- paste0(this_dir, eval_metric,"result/")

library(parallel)
mclapply(1:length(kns), function(i){
  summary_KBET_d6(meta_ls_13, meta_ls_23, data_dir, methods_use, method_dir, dataset_use, eval_metric, this_dir, kns[i])
}, mc.cores = 5)



pct = c(seq(from = 5, to = 25, by = 5))
kns = floor(pct *  mean_sample_size / 100)
uc <- paste0(pct, "%")
for(i in 1:length(kns)){
  kn = kns[i]
  dir_result <- paste0(res_dir, kn, "/kbet_median_acceptance_rate_total.csv")

  if(file.exists(dir_result)){
    res_i = read.csv(dir_result, header=T)
    print(colnames(res_i))
    rej_i <- subset(res_i, select = c("methods","median_observed"))
    acc_i <- subset(res_i, select = c("methods","median_acceptance_rate"))
    print(dim(rej_i))
    if(i==1){
      rej = rej_i
      acc = acc_i
    }else{
      rej = merge(rej, rej_i, by="methods")  
      acc = merge(acc, acc_i, by="methods") 
    }
  }  
}

colnames(rej) <- c("methods", uc)
colnames(acc) <- c("methods", uc)
med_rej <- c()
med_acc <- c()
for(r in 1:nrow(rej)){
  med_rej <- c(med_rej, median(as.numeric(rej[r,uc])))
  med_acc <- c(med_acc, median(as.numeric(acc[r,uc])))
}
rej$median_val <- med_rej
acc$median_val <- med_acc

rej = rej[order(rej$median_val, decreasing = F),]
acc = acc[order(acc$median_val, decreasing = T),]

write.csv(rej, paste0(res_dir,dataset_use,"_kbet_rejection.csv"))
write.csv(acc,paste0(res_dir,dataset_use,"_kbet_acception.csv"))


