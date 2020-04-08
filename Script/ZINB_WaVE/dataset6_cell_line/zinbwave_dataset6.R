source('../ZINB_WaVE_analysis.R')

nThread = 8
base_name = basename(getwd())

# #### preprocessing ####
# b1_rc = fast_read_table('b1_exprs.txt')
# # b1_rc = expm1(b1_rc)
# # ez_head(b1_rc)
# b1_cell_anno = fast_read_table('./b1_celltype.txt')
# b2_rc = fast_read_table('b2_exprs.txt')
# # b2_rc = expm1(b2_rc)
# # ez_head(b2_rc)
# b2_cell_anno = fast_read_table('b2_celltype.txt')
# 
# b3_rc = fast_read_table('b3_exprs.txt')
# # b2_rc = expm1(b2_rc)
# # ez_head(b2_rc)
# b3_cell_anno = fast_read_table('b3_celltype.txt')
# 
# all_b_rc = ezcbind(b1_rc, b2_rc, b3_rc)
# # all_b_rc = round(all_b_rc)
# ### filter genes
# # gene_exp = rowSums(all_b_rc)
# # # hist(log1p(gene_exp))
# # all_b_rc = all_b_rc[gene_exp > 0, ]
# # # ez_head(all_b_rc)
# 
# all_b_cell_anno = ezrbind(b1_cell_anno, b2_cell_anno, b3_cell_anno)
# temp =  create_group(list(colnames(b1_rc), colnames(b2_rc), colnames(b3_rc)), c("batch1", "batch2", "batch3"))
# all_b_cell_anno$batch = temp[, 1]
# colnames(all_b_cell_anno)[1] = "cell_type"
# 
# all_b_rc_f = seurat_filter_norm_data(all_b_rc)
# all_b_rc = all_b_rc[rowname(all_b_rc_f), colnames(all_b_rc_f), drop = F]
# all_b_cell_anno = all_b_cell_anno[colnames(all_b_rc_f), , drop = F]
# 
# all_b_rc_f = seurat_filter_norm_data(all_b_rc)
# all_b_rc = all_b_rc[rownames(all_b_rc_f), colnames(all_b_rc_f), drop = F]
# all_b_cell_anno = all_b_cell_anno[colnames(all_b_rc_f), , drop = F]
# 
# saveRDS(all_b_rc, "all_filtered_hvg_cells_rc.RDS")
# saveRDS(all_b_cell_anno, "all_filtered_hvg_cells_cell_anno.RDS")
# #######################

all_b_rc = fast_read_table('./all_filtered_hvg_cells_rc.RDS')
all_b_cell_anno = fast_read_table('./all_filtered_hvg_cells_cell_anno.RDS')

ZINB_WaVE_analysis(all_b_rc, all_b_cell_anno, base_name = base_name, nthread = nThread, bigdata = T, prop_fit = 0.1)



