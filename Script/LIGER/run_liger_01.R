
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
outfile_prefix = "Dataset1"
save_obj = F

src_dir = "./"
working_dir = "../../Output/"
read_dir = "../../Data/dataset1/"

expr_filename = 'dataset1_sm_uc3.txt'
metadata_filename = 'sample_sm_uc3.txt'

batch_label = "batchlb"
celltype_label = "CellType"

########################
# read data 

expr_mat <- read.table(file = paste0(read_dir,expr_filename),sep="\t",header=T,row.names=1,check.names = F)
metadata <- read.table(file = paste0(read_dir,metadata_filename),sep="\t",header=T,row.names=1,check.names = F)

colnames(metadata)[colnames(metadata) == 'batch'] <- 'batchlb'
metadata$batch <- ifelse(metadata$batchlb == 'Batch1', 1, 2 )

colnames(metadata)[colnames(metadata) == 'celltype'] <- 'CellType'

expr_mat <- expr_mat[, rownames(metadata)]

########################
# run pipeline

source(paste0(src_dir,'call_liger.R'))

liger_obj <- liger_preprocess(expr_mat, metadata, 
                              var.thresh=var.thresh,
                              batch_label = batch_label)
  
call_liger(liger_obj, metadata, batch_label, celltype_label, k = k, nrep = nrep, 
           plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)


