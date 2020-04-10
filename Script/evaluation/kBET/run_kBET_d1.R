library(ggplot2)
library(kBET)
library(ggpubr)
library(ezTools)

rm(list=ls())
this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset1_uc3/'
setwd(this_dir)
eval_metric <- 'kBET_Hoa_v2/'
dataset_use <- 'dataset1'

utils_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/source_code_github/evaluation_utils/'
source(paste0(utils_dir,'evaluation_utils.R'))
source(paste0(utils_dir,'kbet/kbet_utils.R'))

data_dir = '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/dataset1/'

method_dir <- c('raw','seurat2','seurat3',
                'harmony','fastMNN','classicMNN',
                'combat','limma','scgen',
                'scanorama','resnet','zinbwave',
                'scmerge','liger')

method_dir <- c('scgen')
methods_use <- c('scGen')
methods_use <- c('Raw','Seurat_2','Seurat_3',
                 'Harmony','fastMNN','MNN_Correct',
                 'ComBat','limma','scGen',
                 'Scanorama','MMD-ResNet','ZINB-WaVE',
                 'scMerge','LIGER')


# fn <- 'dataset2_resnet_pca.csv'
# myResnet <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
# head(myResnet)
# rownames(myResnet) <- gsub('-[0-9]$','',rownames(myResnet))
# write.csv(myResnet, file = paste0(data_dir, fn), row.names = T, quote = F)


dir.create(paste0(this_dir, eval_metric))
fn <- paste0(dataset_use, '_raw_pca.csv')
meta_ls <- get_celltype_common_kbet_V2(data_dir, fn)

size = length(meta_ls$cells_extract)
print('******Sample size is: ********')
print(size) 
pct = c(seq(from = 5, to = 25, by = 5))
kns = floor(pct*size/100)
print(kns)

res_dir = paste0(this_dir, eval_metric,'result/')
dir.create(res_dir, showWarnings = F)

for(i in 1:length(kns)){
  
  kn = kns[i]
  summary_KBET_v2(meta_ls, data_dir, methods_use, method_dir, dataset_use, eval_metric, this_dir, kn)
}

generate_median_kbet(pct, kns, res_dir, dataset_use, this_dir, kbet_output=T)

