
# Author : Kok Siong Ang 
# Date : 20/09/2019
# Proj : Run Seurat 3 pipeline

########################
#load packages

#library(Seurat)  # Seurat 2 version
#detach("package:Seurat", unload=TRUE)
library(Seurat, lib.loc = "/scbio4/tools/R/R-3.6.0/library") #R 3.6.0, Seurat v3
packageVersion('Seurat')
library(magrittr)
library(cowplot)

rm(list=ls())

########################
#settings

normData = T
Datascaling = T
regressUMI = F
min_cells = 0
min_genes = 0
norm_method = "LogNormalize"
scale_factor = 10000
numVG = 300
nhvg = 2000
npcs = 20
visualize = F
outfile_prefix = "Dataset8"
save_obj = F

src_dir = "./"
working_dir = "./"
#working_dir = "../../Output/"
read_dir = "/batch_norm/datasets/dataset8_Mouse_brain/filtered_genes_and_cells/"
#read_dir = "../../Data/"

expr_filename = "dropviz_and_nuclei_combined_filtered_UMI.RDS";
cellinfo_filename = "dropviz_and_nuclei_combined_filtered_cell_info.txt"

batch_label = "batchlb"
celltype_label = "CellType"

########################
# read data 

expr_mat = readRDS(paste0(read_dir,expr_filename))
metadata = read.table(file = paste0(read_dir,cellinfo_filename),sep="\t",header=T,row.names=1,check.names = F)

expr_mat <- expr_mat[, rownames(metadata)]

colnames(metadata)[colnames(metadata) == 'batch'] <- batch_label
metadata$batch <- ifelse(metadata[, batch_label] == 'batch1', 1, 2)
colnames(metadata)[colnames(metadata) == 'cell_type'] <- celltype_label

########################
# run pipeline

source(paste0(src_dir,'call_seurat_3.R'))
#setwd(working_dir)

batch_list = seurat3_preprocess(
                expr_mat, metadata, 
                normData = normData, Datascaling = Datascaling, regressUMI = regressUMI, 
                min_cells = min_cells, min_genes = min_genes, 
                norm_method = norm_method, scale_factor = scale_factor, 
                numVG = numVG, nhvg = nhvg, 
                batch_label = batch_label, celltype_label = celltype_label)

batches = call_seurat3(batch_list, batch_label, celltype_label, npcs, plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)




