source('../ZINB_WaVE_analysis.R')
b1_rc = fast_read_table('b1c_exprs.txt')
# b1_rc = expm1(b1_rc)
# ez_head(b1_rc)
b1_cell_anno = fast_read_table('./b1c_celltype.txt')
b2_rc = fast_read_table('b2c_exprs.txt')
# b2_rc = expm1(b2_rc)
# ez_head(b2_rc)
b2_cell_anno = fast_read_table('b2c_celltype.txt')

all_b_rc = ezcbind(b1_rc, b2_rc)
# all_b_rc = round(all_b_rc)
### filter genes
gene_exp = rowSums(all_b_rc)
# hist(log1p(gene_exp))
all_b_rc = all_b_rc[gene_exp > 0, ]
# ez_head(all_b_rc)

all_b_cell_anno = ezrbind(b1_cell_anno, b2_cell_anno)
temp =  create_group(list(colnames(b1_rc), colnames(b2_rc)), c("batch1", "batch2"))
all_b_cell_anno$batch = temp[, 1]
colnames(all_b_cell_anno)[1] = "cell_type"

ZINB_WaVE_analysis(all_b_rc, all_b_cell_anno, "Dataset7_mouse_retina")



