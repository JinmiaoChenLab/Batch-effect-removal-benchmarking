
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
outfile_prefix = "Dataset6"
save_obj = F

src_dir = "./"
working_dir = "../../Output/"
read_dir = "../../Data/dataset6/"

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

liger_obj <- liger_preprocess(expr_mat, metadata, 
                              var.thresh=var.thresh,
                              batch_label = batch_label)

call_liger(liger_obj, metadata, batch_label, celltype_label, k = k, nrep = nrep, 
           plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)
