source('../ZINB_WaVE_analysis.R')
nThread = 8
base_name = basename(getwd())

# #### preprocessing ####
# b1_rc = fast_read_table('dataset10_Mouse_BM_hematopoietic/b1_exprs.txt')
# b1_cell_anno = fast_read_table('dataset10_Mouse_BM_hematopoietic/b1_celltype.txt')
# b2_rc = fast_read_table('dataset10_Mouse_BM_hematopoietic/b2_exprs.txt')
# b2_cell_anno = fast_read_table('dataset10_Mouse_BM_hematopoietic/b2_celltype.txt')
# 
# all_b_rc = ezcbind(b1_rc, b2_rc)
# all_b_cell_anno = ezrbind(b1_cell_anno, b2_cell_anno)
# temp =  create_group(list(colnames(b1_rc), colnames(b2_rc)), c("batch1", "batch2"))
# all_b_cell_anno$batch = temp[, 1]
# 
# a = fast_read_table('../../dataset/dataset10_hematoMNN_Hoa/dataset10_hematoMNN_Hoa_all_batch_rc.txt')
# 
# #######################

all_b_rc = fast_read_table('../../dataset/dataset10_hematoMNN_Hoa/dataset10_hematoMNN_Hoa_all_batch_rc.txt')
all_b_cell_anno = fast_read_table('../../dataset/dataset10_hematoMNN_Hoa/dataset10_hematoMNN_Hoa_all_cell_anno.txt')

ZINB_WaVE_analysis(all_b_rc, all_b_cell_anno, base_name = base_name, nthread = nThread, bigdata = F)

