# Run all kBET

library(ggplot2)
library(kBET)
library(ggpubr)
library(ezTools)


##################
## dataset 10
##################

rm(list=ls())
this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset10_hematoMNN_Hoa/'
setwd(this_dir)
data_dir = '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/dataset10/'
eval_metric <- 'kBET_Hoa_v2/'
dataset_use <- 'dataset10'

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
size = length(meta_ls$cells_extract)
print('******Sample size is: ********')
print(size)  # 1626
pct = c(seq(from = 5, to = 25, by = 5))
kns = floor(pct * size / 100)

res_dir <- paste0(this_dir, eval_metric,"result/")



library(parallel)
mclapply(1:length(kns), function(i){
  kn = kns[i]
  print(paste0('******* k0 value use is: ', kn))
  summary_KBET_v2(meta_ls, data_dir, methods_use, method_dir, dataset_use, eval_metric, this_dir, kn)
  
}, mc.cores = 4)


generate_median_kbet(pct, kns, res_dir, dataset_use, this_dir, kbet_output=T)
