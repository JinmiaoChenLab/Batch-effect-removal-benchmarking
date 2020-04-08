source('../ZINB_WaVE_analysis.R')

base_name = basename(getwd())
load('scaled_vargens_filtered_batch1_seqwell_batch2_10x.rda')
# b1_rc = expm1(b1_rc)
# ez_head(b1_rc)
# b1_cell_anno = read.table('./sample_sm_uc3.txt')
# b2_rc = fast_read_table('b2_exprs.txt')
# # b2_rc = expm1(b2_rc)
# # ez_head(b2_rc)
# b2_cell_anno = read.table('b2_celltype.txt')

all_b_rc = pbmc@raw.data
# all_b_rc = round(all_b_rc)
### filter genes
gene_exp = rowSums(all_b_rc)
# hist(log1p(gene_exp))
all_b_rc = all_b_rc[gene_exp > 0, ]
# ez_head(all_b_rc)

all_b_cell_anno = pbmc@meta.data
# temp =  create_group(list(colnames(b1_rc), colnames(b2_rc)), c("batch1", "batch2"))
# all_b_cell_anno$batch = temp[, 1]
colnames(all_b_cell_anno)[colnames(all_b_cell_anno) == "ct"] = "cell_type"

ZINB_WaVE_analysis(all_b_rc, all_b_cell_anno, base_name)



