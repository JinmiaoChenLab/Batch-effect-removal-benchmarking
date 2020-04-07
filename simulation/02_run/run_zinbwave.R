
this.dir <- '/home/marion/Marion/Project/Hoa_batch_normalization/simulation_dataset_V3/'
setwd(this.dir)

source('demo_zinb_wave/ZINB_WaVE_analysis.R')
library(fs)
library(parallel)
library(ezTools)
library(cowplot)

nThread = 25

#' This script is to remove batcheffect using ZINB_WaVE
#' 
#' Author:Xiaomeng
#' Date: Fri Jul  5 08:35:20 2019 ------------------------------
library(zinbwave)
library(scRNAseq)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)
library(tictoc)
library(ezTools)
library(cowplot)
library(lubridate)
library(dplyr)

ZINB_WaVE_analysis = function(expression_data, cell_anno, base_name = 'project', K = 50, nthread = 8
                              , bigdata = F, prop_fit = 0.1){
  ez_print("Start ", base_name, " zinbwave analysis")
  ez_print("Data dim: ", dim(expression_data))
  ez_print("Num thread ", nthread)
  if(bigdata){
    ez_print("Using zinbsurf. Prop_fit = ", prop_fit)
  }
  tic(paste0("ZINB_WaVE_", base_name))
  
  colData <- DataFrame(batch = cell_anno[, "batch"],
                       row.names=rownames(cell_anno))
  all_b =   SummarizedExperiment(assays = list(counts = as.matrix(expression_data)), colData = colData)
  
  if(!bigdata){
    b_cov <- zinbwave(all_b, K=K, X="~batch", epsilon=1000, verbose = T, BPPARAM = MulticoreParam(nthread)
                      , normalizedValues=TRUE, residuals = TRUE)
  } else {
    b_cov <- zinbsurf(all_b,  K=K, X="~batch", epsilon=1000, verbose = T, BPPARAM = MulticoreParam(nthread)
                      , prop_fit = prop_fit)
  }
  W <- reducedDim(b_cov)
  
  fast_save_table(W, base_name, "_zinbwave_corrected.txt")
  fast_save_table(b_cov@assays[['normalizedValues']], base_name, "_zinbwave_corrected_matrix.txt")
  saveRDS(b_cov, paste0(base_name, "_after_zinbwave_corrected.RDS"))
  
  pt_size = 0.5
  if(nrow(expression_data) > 10000){
    pt_size = 0.1
  }
  set_ggplot_settings("pointsize", pt_size)
  set_ggplot_settings("page_width", 10)
  set_ggplot_settings("page_height", 4)
  p1 = plot_scatter(W, color = cell_anno[, "batch", drop = F], title = paste0(base_name,  " After ZINB_WaVE"))
  p2 = plot_scatter(W, color = cell_anno[, "cell_type", drop = F], title = paste0(base_name,  " After ZINB_WaVE"))
  p = plot_grid(p1, p2)
  save_image(p, paste0(base_name, " After ZINB_WaVE.png"))
  
  b_cov <- scater::runTSNE(b_cov, exprs_values = "normalizedValues")
  colnames(b_cov@reducedDims$TSNE) <- c('tSNE_1','tSNE_2')
  rownames(b_cov@reducedDims$TSNE) <- colnames(b_cov@assays[["normalizedValues"]])
  tsne_df <- merge(b_cov@reducedDims$TSNE,cell_anno,by=0)
  tsne_df <- tsne_df[,c('tSNE_1','tSNE_2','batch','cell_type')]
  fast_save_table(tsne_df, base_name, "_tsne.txt")
  
  p1 = plot_scatter(tsne_df, color = tsne_df$batch, title = paste0(base_name,  " After ZINB_WaVE"))
  p2 = plot_scatter(tsne_df, color = tsne_df$cell_type, title = paste0(base_name,  " After ZINB_WaVE"))
  p = plot_grid(p1, p2)
  save_image(p, paste0(base_name, " After ZINB_WaVE_from_norm.png"))

  return(NULL)
}


# #### preprocessing ####
base_name = basename(getwd())

files = dir('data', '.*0$', full.names = T)

all_rc_and_anno_list = lapply(1:length(files), function(i){
  if(is_file(files[i])){
    return(NULL)
  }
  print(files[i])
  data_dir = files[i]
  # geneinfo <- read.table(paste0(data_dir,"/geneinfo.txt"), head=T, row.names = 1, check.names = FALSE)
  cellinfo <- read.table(paste0(data_dir,"/cellinfo.txt"), head=T, row.names = 1, check.names = FALSE)
  counts <- read.table(paste0(data_dir,"/counts.txt"), head=T, row.names = 1, check.names = FALSE)

  counts_hm <- t(counts)
  dim(counts_hm)
  ez_head(counts_hm)
  colnames(counts_hm) <- paste(cellinfo[colnames(counts_hm),'Group'],colnames(counts_hm),sep="_")
  rownames(cellinfo) <- paste(cellinfo[rownames(cellinfo),'Group'],rownames(cellinfo), sep="_")
  colnames(cellinfo)[2:3] = c("batch", "cell_type")
  ez_head(cellinfo)
  cellinfo <- cellinfo[colnames(counts_hm), , drop = F]
  res = list(counts = counts_hm, cell_anno = cellinfo)
  return(res)
})
names(all_rc_and_anno_list) = basename(files)
all_rc_and_anno_list2 = all_rc_and_anno_list[!sapply(all_rc_and_anno_list, is.null)]
saveRDS(all_rc_and_anno_list2, "all_rc_and_anno_list.RDS")

all_rc_and_anno_list = fast_read_table('/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/ZINB_WaVE/dataset3_simulation_v2/analysis/all_rc_and_anno_list.RDS')

lapply(1:length(all_rc_and_anno_list), function(i){
  print(names(all_rc_and_anno_list)[i])
  ZINB_WaVE_analysis(all_rc_and_anno_list[[i]]$counts
                     , all_rc_and_anno_list[[i]]$cell_anno
                     , base_name = names(all_rc_and_anno_list)[i]
                     , nthread = nThread)
})

HVG_rc_and_anno_list = fast_read_table('demo_zinb_wave/HVG_rc_and_anno_list.RDS')

lapply(1:length(HVG_rc_and_anno_list), function(i){
  print(names(all_rc_and_anno_list)[i])
  ZINB_WaVE_analysis(all_rc_and_anno_list[[i]]$counts
                     , all_rc_and_anno_list[[i]]$cell_anno
                     , base_name = names(all_rc_and_anno_list)[i]
                     , nthread = nThread)
})


