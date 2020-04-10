
library(ggplot2)
library(kBET)
library(ggpubr)
library(ezTools)


##################
## dataset 9
##################

rm(list=ls())
this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset9_Human_cell_atlas/'
setwd(this_dir)
data_dir = '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/dataset9/'
# eval_metric <- 'kBET_Hoa_ds_v2/' 
eval_metric <- 'kBET_Hoa_ds_02/' 
# kn = 30
dataset_use <- 'dataset9'

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
percent_ext <- 0.2
meta_ls <- get_downsample_d9(data_dir, fn, percent_ext)


#### iteration 
size = length(meta_ls$cells_extract)  # should be 6000, 10 % cells


# for our evaluation
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
  summary_KBET_d9(meta_ls, data_dir, methods_use, method_dir, dataset_use, eval_metric, this_dir, kn)
}, mc.cores = 8)


###################
## kbet generating results
###################

generate_median_kbet(pct, kns, res_dir, dataset_use, this_dir, kbet_output=T)


