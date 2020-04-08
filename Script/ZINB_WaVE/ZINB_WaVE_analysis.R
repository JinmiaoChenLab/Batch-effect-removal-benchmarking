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
  
  # b_ori = zinbwave(all_b, K=2, epsilon=1000, verbose = T)
  # W2 <- reducedDim(b_ori)
  # 
  # p = plot_scatter(W2, color = cell_anno[, "batch", drop = F], title =  paste0(base_name,  " Before ZINB_WaVE"))
  # save_image(p, "Dataset10 Before ZINB_WaVE.png")
  # 
  # save(b1_rc, b1_cell_anno, b2_rc, b2_cell_anno, expression_data, cell_anno, colData, all_b, b_cov, W, b_ori, W2,
  #      file = paste0("ZIB_WaVE_", base_name, ".Rdata"))
  # return(list(befor_zb = b_ori, after_zb = b_cov))
  
  b = toc()
  p = (b$toc - b$tic) %>% seconds_to_period() %>% as.character()
  time = init_sample_info(c("Running time", "dim", "thread_num", "Use zinbsurf", "prop_fit")
                          , c(p, paste(dim(expression_data), collapse = "x"), nthread
                              , bigdata, prop_fit), name = b$msg)
  fast_save_table(time, base_name, "_ZINB_WaVE_running_time.txt")
  return(NULL)
}











