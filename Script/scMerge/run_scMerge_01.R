
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
min_cells = 5
min_genes = 300
cosineNorm = FALSE
kmeansK = c(3,3)
replicate_prop = 0.5
npcs=20
visualize = T
outfile_prefix = "Dataset1"
save_obj = F

data("segList", package = "scMerge") 
seg = segList$human$human_scSEG

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

source(paste0(src_dir,'call_scMerge.R'))

scMerge_obj <- scMerge_preprocess(expr_mat, metadata, 
                                  batch_label = batch_label,
                                  filter_genes = filter_genes, filter_cells = filter_cells,
                                  min_cells = min_cells, min_genes = min_genes, 
                                  cosineNorm = cosineNorm)

call_scMerge(scMerge_obj, batch_label, celltype_label, seg = seg, 
             kmeansK = kmeansK, replicate_prop=replicate_prop, npcs=npcs,
             plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)
