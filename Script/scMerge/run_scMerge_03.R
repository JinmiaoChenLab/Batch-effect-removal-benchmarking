
# Author : Marion Chevrier
# Date : 15/10/2019
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
kmeansK = c(2,2)
replicate_prop = 0.5
npcs=20
visualize = T
outfile_prefix = "Dataset3"
save_obj = F

src_dir = "./"
read_dir = "../../Data/dataset3/"

lsdir <- paste0(list.dirs(read_dir, recursive=FALSE),'/')

sapply(lsdir,function(x){
  
  x2 <- gsub(paste0(read_dir,'/'),'',x)
  
  selection <- c('HVG','all')
  
  sapply(selection, function(s){
    
    if(s=='HVG'){
      expr_mat_filename = 'counts_HVG.txt'
    } else {
      expr_mat_filename = 'counts.txt'
    }
    
    metadata_filename = 'cellinfo.txt'
    
    batch_label = "batchlb"
    celltype_label = "CellType"
    
    working_dir = paste0("../../Output/Dataset3/",x2,s,"/")
    dir.create("../../Output/Dataset3/", showWarnings = FALSE)
    dir.create(paste0("../../Output/Dataset3/",x2), showWarnings = FALSE)
    dir.create(paste0("../../Output/Dataset3/",x2,s,"/"), showWarnings = FALSE)
    
    ########################
    # read data 
    
    expr_mat <- read.table(file = paste0(x,expr_mat_filename),sep="\t",header=T,row.names=1,check.names = F)
    expr_mat <- t(expr_mat)
    metadata <- read.table(file = paste0(x,metadata_filename),sep="\t",header=T,row.names=1,check.names = F)
    
    colnames(metadata)[colnames(metadata) == 'Batch'] <- 'batchlb'
    metadata$batch <- ifelse(metadata$batchlb == 'Batch1', 1, 2 )
    
    colnames(metadata)[colnames(metadata) == 'Group'] <- 'CellType'
    
    expr_mat <- expr_mat[, rownames(metadata)]

    ########################
    # seg
    
    geneinfo <- read.table(paste0(x,'/geneinfo.txt'), head=T, sep='\t')
    geneinfo_wobatch <- geneinfo[geneinfo$DEFacGroup1==1 & geneinfo$DEFacGroup2==1,]
    counts2 <- expr_mat[rownames(expr_mat) %in% as.character(geneinfo_wobatch$Gene),]
    gene_lowvar <- names(sort(apply(counts2,1,var), decreasing=F)[1:200])
    
    ########################
    # run pipeline
    
    source(paste0(src_dir,'call_scMerge.R'))
    
    scMerge_obj <- scMerge_preprocess(expr_mat, metadata, 
                                      batch_label = batch_label,
                                      filter_genes = filter_genes, filter_cells = filter_cells,
                                      min_cells = min_cells, min_genes = min_genes, 
                                      cosineNorm = cosineNorm)
    
    call_scMerge(scMerge_obj, batch_label, celltype_label, seg = gene_lowvar,
                 kmeansK = kmeansK, marker = rownames(expr_mat), replicate_prop=replicate_prop, npcs=npcs,
                 plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)
    
  })
  
})
