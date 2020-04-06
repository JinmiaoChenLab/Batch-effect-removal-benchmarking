
# Author : Marion Chevrier
# Date : 30/09/2019
# Proj : Run LIGER pipeline

########################
#load packages

#devtools::install_github('MacoskoLab/liger')
library(scales)
library(liger)
library(Matrix)
library(Rtsne)
library(ggplot2)

rm(list=ls())

########################
#settings

var.thresh = 0.1
k = 20
nrep = 3
visualize = T
outfile_prefix = "Dataset9"
save_obj = F

src_dir = "./"
working_dir = "../../Output/"
read_dir = "../../Data/dataset9/"

expr_filename = "HCA_genes_cells_filtered_filtered_UMI.RDS";
cellinfo_filename = "HCA_genes_cells_filtered_filtered_cell_info.txt"

batch_label = "batchlb"
celltype_label = "CellType"

########################
# read data 

expr_mat = readRDS(paste0(read_dir,expr_filename))
metadata = read.table(file = paste0(read_dir,cellinfo_filename),sep="\t",header=T,row.names=1,check.names = F)

expr_mat <- expr_mat[, rownames(metadata)]

colnames(metadata)[colnames(metadata) == 'tissue'] <- batch_label
metadata$batch <- ifelse(metadata[,batch_label] == 'Cord blood', 1, 2)
colnames(metadata)[colnames(metadata) == 'batch'] <- celltype_label

########################
# run pipeline

source(paste0(src_dir,'call_liger.R'))

liger_obj <- liger_preprocess(expr_mat, metadata, 
                              var.thresh=var.thresh,
                              batch_label = batch_label)

call_liger(liger_obj, metadata, batch_label, celltype_label, k = k, nrep = nrep, 
           plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)
