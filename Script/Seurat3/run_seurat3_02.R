
# Author : Kok Siong Ang 
# Date : 03/09/2019
# Proj : Run Seurat 3 pipeline

########################
#load packages

#library(Seurat)  # Seurat 2 version
#detach("package:Seurat", unload=TRUE)
library("Seurat", lib.loc="/tools/R/R-3.5.0/library") #Seurat v3
packageVersion('Seurat')
library(magrittr)
library(cowplot)

rm(list=ls())

########################
#settings

normData = T
Datascaling = T
regressUMI = F
min_cells = 0 #10
min_genes = 0 #300
norm_method = "LogNormalize"
scale_factor = 10000
numVG = 300
nhvg = 2000
npcs = 20
visualize = T
outfile_prefix = "Dataset2"
save_obj = F

src_dir = "./"
working_dir = "../../Output/"
read_dir = "../../Data/dataset2/"

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




