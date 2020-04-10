# demo COMBAT 
# author: created by hoa tran modified by jinmiao chen
# Documentation: https://rdrr.io/bioc/sva/man/ComBat.html

rm(list=ls())

source("/home/chenjm/Hoa_batch_normalization/demo_combat/combat_functions.R")

library(sva)
library(Seurat)
library(matrixStats)

packageVersion('Seurat')

this.dir <- '/home/chenjm/Hoa_batch_normalization/demo_combat'
setwd(this.dir)

data_dir <- '/home/xm/Projects/batch_norm/datasets/dataset9_Human_cell_atlas/filtered_genes_and_cells/'

base_name <- 'combat_dataset9/'
dir.create(base_name, showWarnings = FALSE)

TPM_file <- 'HCA_genes_cells_filtered_filtered_UMI.RDS'  # replace by link to dataset
sample_file <- 'HCA_genes_cells_filtered_filtered_cell_info.txt' # replace by link to dataset

myData <- readRDS(paste0(data_dir,TPM_file))
mySample <- read.table(paste0(data_dir,sample_file),sep="\t",header=T,row.names=1,check.names=F)

mySample$celltype <- mySample$cell_type

mySample$batch <- ifelse(mySample$tissue=="Bone marrow","batch1",
                         ifelse(mySample$tissue=="Cord blood","batch2",NA))


print(unique(mySample$batch))

# Pre-filtering data 
# First option: use Seurat workflow to normalize and scale data first, before use as input to Combat function
# myData <- filter_data_seurat(myData, mySample, min_cells=10, min_genes=300, group_col='celltype', regressUMI=FALSE)

# Second option: use filtering function as below
# Func as.matrix() can not work with big data, divide big matrix into small parts and do filtering 

print("Median normalizing counts and log-transforming")
col_sums = Matrix::colSums(myData)
med_trans = median(col_sums)
norm_data = med_trans * myData/col_sums
rm(myData)

myFilteredData <- log(preprocess_big(norm_data)+1,2)
rm(norm_data)

row_max = matrixStats::rowMaxs(myFilteredData)
row_min = matrixStats::rowMins(myFilteredData)

if(sum(row_max==row_min)>0){
  myFilteredData = myFilteredData[row_max != row_min,]
}

cells_use <- colnames(myFilteredData)
mySample <- mySample[cells_use,]

# save filtered data for combat and limma
#saveRDS(myFilteredData, file = paste0(base_name,"combat_filtered_mtx.rds"))
#saveRDS(mySample, file = paste0(base_name,"mySample_filtered_cells.rds"))

#########################################
### Remove no variance genes
##
#########################################
# For COMBAT, if applying this method for matrix in which genes with no variance, 
# the function will display error
# So we remove the genes with no variance first, before apply COMBAT

# Check genes variance and remove genes with no variance
# length(which ((apply(pbmc@raw.data, 1, var)!=0) == "FALSE"))
# length(which ((apply(pbmc@raw.data, 1, var)!=0) == "TRUE"))
# rv_genes <- which ((apply(pbmc@raw.data, 1, var)!=0) == "FALSE")
# rv_genes_names <- rownames(pbmc@raw.data)[rv_genes]

#rv_genes <- which ((apply(myFilteredData, 1, var)!=0) == "FALSE")   
#rv_genes_names <- rownames(myFilteredData)[rv_genes]             

# remove no variance genes
#if(rv_genes_names) {
#  myFilteredData <- myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
#} else{
#  print('All genes have variances, continue to run Combat')
#}
#dim(myFilteredData)



#########################################
### COMBAT
#########################################


# First option, opt1 <- 'parametric/'
# mod=NULL, par.prior=TRUE, prior.plots=FALSE
opt1 <- './combat_dataset9/parametric/'
opt2 <- './combat_dataset9/non-parametric/'
base_name <- opt1  #save output in this folder

combat_output <- combat_batch_effect_removal(myFilteredData, batch=mySample$batch, 
                                             mod=NULL, par_prior=TRUE, prior_plots=FALSE, mean_only=FALSE, 
                                             base_name=opt1, save_txt=FALSE, saveRDS=TRUE)

rm(list=ls())

# base_name <- opt2  #save output in this folder
# Second option, opt2 <- 'non-parametric/'
# mod=NULL, par.prior=FALSE,prior.plots=FALSE, mean.only=TRUE
# combat_output <- combat_batch_effect_removal(pbmc@scale.data, batch=pbmc@meta.data[cells_use,'batch'],
#                             mod=NULL, par_prior=FALSE, prior_plots=FALSE, mean_only=TRUE,
#                             base_name=opt2, save_txt=FALSE, saveRDS=TRUE)






