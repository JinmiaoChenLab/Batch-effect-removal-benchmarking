#' This script is to select the cell from koksiong output
#' 
#' Author:Xiaomeng
#' Date: Thu Aug 29 20:35:17 2019 ------------------------------
library(ezTools)

all_rc = read.table('./dataset1_sm_uc3_all.txt', sep = "\t")
sample = read.table('./sample_sm_uc3_all.txt', sep = "\t")
umap = fast_read_table('./harmony_umap.txt')
selected_cells = rownames(umap)

selected_rc = all_rc[, selected_cells, drop = F]
selected_cell_anno = sample[selected_cells, , drop = F]

fast_save_table(selected_rc, "", "dataset1_sm_uc3.txt")
fast_save_table(selected_cell_anno, "", "sample_sm_uc3.txt")

a = fast_read_table('./dataset1_sm_uc3.txt')





