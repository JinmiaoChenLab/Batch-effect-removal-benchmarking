

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
kmeansK = c(1,1,2)
replicate_prop = 0.5
npcs=20
visualize = T
outfile_prefix = "Dataset6"
save_obj = F

data("segList", package = "scMerge") 
seg <- segList$human$human_scSEG

src_dir = "./"
working_dir = "../../Output/"
read_dir = "../../Data/dataset6"

b1_exprs_filename = "b1_exprs.txt"
b2_exprs_filename = "b2_exprs.txt"
b3_exprs_filename = "b3_exprs.txt"
b1_celltype_filename = "b1_celltype.txt"
b2_celltype_filename = "b2_celltype.txt"
b3_celltype_filename = "b3_celltype.txt"

batch_label = "batchlb"
celltype_label = "CellType"

########################
# read data from flat text

b1_exprs <- read.table(file = paste0(read_dir,b1_exprs_filename),sep="\t",header=T,row.names=1,check.names = F)
b2_exprs <- read.table(file = paste0(read_dir,b2_exprs_filename),sep="\t",header=T,row.names=1,check.names = F)
b3_exprs <- read.table(file = paste0(read_dir,b3_exprs_filename),sep="\t",header=T,row.names=1,check.names = F)
b1_celltype <- read.table(file = paste0(read_dir,b1_celltype_filename),sep="\t",header=T,row.names=1,check.names = F)
b2_celltype <- read.table(file = paste0(read_dir,b2_celltype_filename),sep="\t",header=T,row.names=1,check.names = F)
b3_celltype <- read.table(file = paste0(read_dir,b3_celltype_filename),sep="\t",header=T,row.names=1,check.names = F)

b1_celltype$cell <- rownames(b1_celltype)
b1_celltype <- b1_celltype[colnames(b1_exprs),]
b2_celltype$cell <- rownames(b2_celltype)
b2_celltype <- b2_celltype[colnames(b2_exprs),]
b3_celltype$cell <- rownames(b3_celltype)
b3_celltype <- b3_celltype[colnames(b3_exprs),]
b1_metadata <- as.data.frame(b1_celltype)
b2_metadata <- as.data.frame(b2_celltype)
b3_metadata <- as.data.frame(b3_celltype)
b1_metadata$batch <- 1
b2_metadata$batch <- 2
b3_metadata$batch <- 3
b1_metadata$batchlb <- 'Batch_1'
b2_metadata$batchlb <- 'Batch_2'
b3_metadata$batchlb <- 'Batch_3'

expr_mat = cbind(b1_exprs, b2_exprs, b3_exprs)
metadata = rbind(b1_metadata, b2_metadata, b3_metadata)

expr_mat <- expr_mat[, rownames(metadata)]

########################
# run pipeline

source(paste0(src_dir,'call_liger.R'))

source(paste0(src_dir,'call_scMerge.R'))

scMerge_obj <- scMerge_preprocess(expr_mat, metadata, 
                                  batch_label = batch_label,
                                  filter_genes = filter_genes, filter_cells = filter_cells,
                                  min_cells = min_cells, min_genes = min_genes, 
                                  cosineNorm = cosineNorm)

call_scMerge(scMerge_obj, batch_label, celltype_label, seg = seg, 
             kmeansK = kmeansK, replicate_prop=replicate_prop, npcs=npcs,
             plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)
