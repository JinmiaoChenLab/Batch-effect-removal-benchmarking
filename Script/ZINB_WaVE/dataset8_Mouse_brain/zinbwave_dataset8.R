source('../ZINB_WaVE_analysis.R')

nThread = 8
base_name = basename(getwd())

# #### preprocessing ####
# umi = fast_read_table('../../dataset/Mouse_brain/filtered_genes_and_cells/dropviz_and_nuclei_combined_filtered_UMI.txt')
# cell_anno = fast_read_table('../../dataset/Mouse_brain/filtered_genes_and_cells/dropviz_and_nuclei_combined_filtered_cell_info.txt')
# hvg = fast_read_table('../../dataset/Mouse_brain/dropviz_and_nuclei_variable_genes.txt')
# umi2 = umi[hvg[, 1], , drop = F]
# saveRDS(umi2, "dropviz_and_nuclei_combined_filtered_hvg_UMI.RDS")
# # fast_save_table(umi2, "", "../../dataset/Mouse_brain/filtered_genes_and_cells/dropviz_and_nuclei_combined_filtered_hvg_UMI.txt")
# # fast_save_table(umi2, "", "dropviz_and_nuclei_combined_filtered_hvg_UMI.txt")
# #######################



all_b_rc = fast_read_table('../../dataset/Mouse_brain/dropviz_and_nuclei_combined_hvg_filtered_cells_UMI.RDS')
all_b_cell_anno = fast_read_table('../../dataset/Mouse_brain/dropviz_and_nuclei_combined_filtered_cell_info.txt')

ZINB_WaVE_analysis(all_b_rc, all_b_cell_anno, base_name = base_name, nthread = nThread, bigdata = F)
# ZINB_WaVE_analysis(all_b_rc, all_b_cell_anno, base_name = base_name, nthread = nThread, bigdata = T, prop_fit = 0.005)



