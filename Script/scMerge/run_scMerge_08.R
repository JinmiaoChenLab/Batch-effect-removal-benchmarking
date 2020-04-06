
# Author : Marion Chevrier
# Date : 01/10/2019
# Proj : Run scMerge pipeline

########################
#load packages

#BiocManager::install("scMerge")
library(scMerge)
library(SingleCellExperiment)
library(scater)
library(scran)
library(cowplot)
library(Rtsne)
library(stats)

rm(list=ls())

########################
#settings

filter_genes = F
filter_cells = F
cosineNorm = FALSE
kmeansK = c(12,11)
replicate_prop = 0.5
npcs=20
visualize = T
outfile_prefix = "Dataset8"
save_obj = F

data("segList", package = "scMerge") 
seg <- segList$mouse$mouse_scSEG

src_dir = "./"
working_dir = "../../Output/"
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

source(paste0(src_dir,'call_scMerge.R'))

scMerge_obj <- scMerge_preprocess(expr_mat, metadata, 
                                  batch_label = batch_label,
                                  filter_genes = filter_genes, filter_cells = filter_cells,
                                  min_cells = min_cells, min_genes = min_genes, 
                                  cosineNorm = cosineNorm)

call_scMerge(scMerge_obj, batch_label, celltype_label, seg = seg,
             kmeansK = kmeansK, replicate_prop=replicate_prop, npcs=npcs,
             plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)
