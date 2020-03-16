
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
min_cells = 0
min_genes = 0
norm_method = "LogNormalize"
scale_factor = 10000
numVG = 300
nhvg = 2000
npcs = 20
visualize = T
outfile_prefix = "Dataset5"
save_obj = F

src_dir = "./"
working_dir = "../../Output/"
read_dir = "../../Data/dataset5/"

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




