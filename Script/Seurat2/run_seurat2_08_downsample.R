
# Author : Kok Siong Ang 
# Date : 03/09/2019
# Proj : Run Seurat 2 pipeline

########################
#load packages

library(cowplot)
library(Seurat)  # Seurat 2 version
library(magrittr)

rm(list=ls())

########################
#settings

filter_genes = T
filter_cells = F
normData = T
Datascaling = T
regressUMI = F
min_cells = 10
min_genes = 300
norm_method = "LogNormalize"
scale_factor = 10000
b_x_low_cutoff = 0.0125
b_x_high_cutoff = 3
b_y_cutoff = 0.5
numVG = 300
nhvg = 2000
npcs = 20
visualize = T
outfile_prefix = "Dataset8_10"
save_obj = F

src_dir = "./"
working_dir = "../../Output/"
read_dir = "../../Data/dataset8/downsample/"

expr_filename = "downsample.rds";
cellinfo_filename = "downsample_meta.rds"

batch_label = "batchlb"
celltype_label = "CellType"

########################
# read data 

expr_mat = readRDS(paste0(read_dir,expr_filename))
metadata = readRDS(paste0(read_dir,cellinfo_filename))

expr_mat <- expr_mat[, rownames(metadata)]

colnames(metadata)[colnames(metadata) == 'batch'] <- batch_label
metadata$batch <- ifelse(metadata[, batch_label] == 'batch1', 1, 2)
colnames(metadata)[colnames(metadata) == 'cell_type'] <- celltype_label

########################
# run pipeline

source(paste0(src_dir,'call_seurat_2.R'))
#setwd(working_dir)

return_obj = seurat2_preprocess(
                expr_mat, metadata, 
                filter_genes = filter_genes, filter_cells = filter_cells,
                normData = normData, Datascaling = Datascaling, regressUMI = regressUMI, 
                min_cells = min_cells, min_genes = min_genes, 
                norm_method = norm_method, scale_factor = scale_factor, 
                b_x_low_cutoff = b_x_low_cutoff, b_x_high_cutoff = b_x_high_cutoff, b_y_cutoff = b_y_cutoff, 
                numVG = numVG, nhvg = nhvg, 
                batch_label = batch_label, celltype_label = celltype_label)
batch_list = return_obj[[1]]
hvg_union = return_obj[[2]]

b_seurat = call_seurat2(batch_list, hvg_union, batch_label, celltype_label, npcs, plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)




