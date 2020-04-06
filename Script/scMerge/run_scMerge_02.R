
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
cosineNorm = TRUE
kmeansK = c(11,11)
replicate_prop = 0.5
npcs=20
visualize = T
outfile_prefix = "Dataset2"
save_obj = F

data("segList", package = "scMerge") 
seg <- segList$mouse$mouse_scSEG

src_dir = "./"
working_dir = "../../Output/"
read_dir = "../../scbio4/home/koksiong/data/02/"
#read_dir = "../../Data/"

expr_filename = 'filtered_total_batch1_seqwell_batch2_10x.txt'
metadata_filename = 'filtered_total_sample_ext_organ_celltype_batch.txt'

batch_label = "batchlb"
celltype_label = "CellType"

########################
# read data 

expr_mat <- read.table(file = paste0(read_dir,expr_filename),sep="\t",header=T,row.names=1,check.names = F)
metadata <- read.table(file = paste0(read_dir,metadata_filename),sep="\t",header=T,row.names=1,check.names = F)

colnames(metadata)[colnames(metadata) == 'ct'] <- 'CellType'

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
             kmeansK = kmeansK, replicate_prop = replicate_prop, npcs = npcs,
             plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)
