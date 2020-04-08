source('../ZINB_WaVE_analysis.R')

nThread = 8
base_name = basename(getwd())

# #### preprocessing ####
# umi = fast_read_table('../../dataset/dataset9_Human_cell_atlas/filtered_genes_and_cells/HCA_genes_cells_filtered_filtered_UMI.txt')
# cell_anno = fast_read_table('../../dataset/dataset9_Human_cell_atlas/HCA_genes_cells_filtered_filtered_cell_info.txt')
# hvg = fast_read_table('../../dataset/dataset9_Human_cell_atlas/HCA_variable_genes.txt')
# umi2 = umi[hvg[, 1], , drop = F]
# saveRDS(umi2, "HCA_genes_cells_filtered_hvg_filtered_UMI.RDS")
# #######################

all_b_rc = fast_read_table('./HCA_genes_cells_filtered_hvg_filtered_UMI.RDS')
all_b_cell_anno = fast_read_table('./../../dataset/dataset9_Human_cell_atlas/HCA_genes_cells_filtered_filtered_cell_info.txt')

ZINB_WaVE_analysis(all_b_rc, all_b_cell_anno, base_name = base_name, nthread = nThread, bigdata = T, prop_fit = 0.005)



