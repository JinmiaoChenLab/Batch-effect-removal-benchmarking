
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

filter_genes = T
filter_cells = T
min_cells = 10
min_genes = 300
cosineNorm = FALSE
kmeansK = c(5,12)
replicate_prop = 0.5
npcs=20
visualize = T
outfile_prefix = "Dataset7"
save_obj = F

data("segList", package = "scMerge") 
seg <- segList$mouse$mouse_scSEG

src_dir = "./"
working_dir = "../../Output/"
read_dir = "../../Data/dataset7"

b1_exprs_filename = "b1_exprs.txt"
b2_exprs_filename = "b2_exprs.txt"
b1_celltype_filename = "b1_celltype.txt"
b2_celltype_filename = "b2_celltype.txt"

batch_label = "batchlb"
celltype_label = "CellType"

########################
# read data 

b1_exprs <- read.table(file = paste0(read_dir,b1_exprs_filename),sep="\t",header=T,row.names=1,check.names = F)
b2_exprs <- read.table(file = paste0(read_dir,b2_exprs_filename),sep="\t",header=T,row.names=1,check.names = F)
b1_celltype <- read.table(file = paste0(read_dir,b1_celltype_filename),sep="\t",header=T,row.names=1,check.names = F)
b2_celltype <- read.table(file = paste0(read_dir,b2_celltype_filename),sep="\t",header=T,row.names=1,check.names = F)

b1_celltype$cell <- rownames(b1_celltype)
b1_celltype <- b1_celltype[colnames(b1_exprs),]
b2_celltype$cell <- rownames(b2_celltype)
b2_celltype <- b2_celltype[colnames(b2_exprs),]
b1_metadata <- as.data.frame(b1_celltype)
b2_metadata <- as.data.frame(b2_celltype)
b1_metadata$batch <- 1
b2_metadata$batch <- 2
b1_metadata$batchlb <- 'Batch_1'
b2_metadata$batchlb <- 'Batch_2'

expr_mat = cbind(b1_exprs,b2_exprs)
metadata = rbind(b1_metadata, b2_metadata)

expr_mat <- expr_mat[, rownames(metadata)]

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
