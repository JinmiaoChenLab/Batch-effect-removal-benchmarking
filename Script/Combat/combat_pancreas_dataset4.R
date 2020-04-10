# demo COMBAT 
# author: hoa tran
# Documentation: https://rdrr.io/bioc/sva/man/ComBat.html
  
rm(list=ls())
library(sva)
library(Seurat)
packageVersion('Seurat')
library(matrixStats)
this.dir <- '/home/hoa/hoatran/demo_normalization/michelle/demo_combat/'
setwd(this.dir)
data_dir <- '/home/hoa/hoatran/demo_normalization/dataset/dataset4_human_pancreatic/final/'
base_name <- 'combat_dataset4/'
dir.create(base_name, showWarnings = FALSE)
  
  
# load evaluation & preprocess utils files #
utils_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/michelle/'
source(paste0(utils_dir,'combat_functions.R'))


TPM_file <- 'myData_pancreatic_5batches.txt'  # replace by link to dataset
sample_file <- 'mySample_pancreatic_5batches.txt' # replace by link to dataset
myData <- read.table(paste0(data_dir,TPM_file),sep="\t",header=T,row.names=1,check.names=F)
mySample <- read.table(paste0(data_dir,sample_file),sep="\t",header=T,row.names=1,check.names=F)
mySample$batch <- mySample$batchlb

opt1 <- './combat_dataset4/parametric/'
# opt2 <- './combat_dataset4/non-parametric/'


# Pre-filtering data 
# Func as.matrix() can not work with big data, divide big matrix into small parts and do filtering
# Normalize each library to the median of the transcript counts across all cells 
# Then, log transform expression values 
myFilteredData <- filter_data_mtx(myData, base_name, is_filter_cells=TRUE,min_genes=300, 
                                  is_filter_genes=TRUE, min_cells=10)

row_max = matrixStats::rowMaxs(myFilteredData)
row_min = matrixStats::rowMins(myFilteredData)

if(sum(row_max==row_min)>0){
  myFilteredData = myFilteredData[row_max != row_min,]
}

cells_use <- colnames(myFilteredData)
mySample <- mySample[cells_use,]

# For COMBAT, if applying this method for matrix in which genes with no variance, 
# the function will display error
# So we remove the genes with no variance first, before apply COMBAT
  
# Check genes variance and remove genes with no variance
# length(which ((apply(pbmc@raw.data, 1, var)!=0) == "FALSE"))
# length(which ((apply(pbmc@raw.data, 1, var)!=0) == "TRUE"))
# rv_genes <- which ((apply(pbmc@raw.data, 1, var)!=0) == "FALSE")
# rv_genes_names <- rownames(pbmc@raw.data)[rv_genes]
  # rv_genes <- which ((apply(pbmc@scale.data, 1, var)!=0) == "FALSE")   
  # rv_genes_names <- rownames(pbmc@scale.data)[rv_genes]             
  # 
  # # remove no variance genes
  # if(rv_genes_names) {
  #   pbmc@scale.data <- pbmc@scale.data[!(rownames(pbmc@scale.data) %in% rv_genes_names),]
  # } else{
  #   print('All genes have variances, continue to run Combat')
  # }
  # dim(pbmc@scale.data)
  
## COMBAT function 
# Opt1: parametric option: mod=NULL, par.prior=TRUE, prior.plots=FALSE, mean_only=FALSE  the best option
# Opt2: non-parametric adjustment, mean-only version mod=NULL, par.prior=FALSE,prior.plots=FALSE, mean.only=TRUE
# Pls load evaluation functions first to have time_execution function
  
# myData: data frame or matrix, genes x cells 
# batch: a vector of bach label, ex: batch <- pbmc@meta.data[cells_use,'batch']
# In case big data, can save as rda or rds format, save_txt=FALSE
  # combat_batch_effect_removal <- function(myData, batch=NULL, mod=NULL, par_prior=TRUE, prior_plots=FALSE, mean_only=FALSE, 
  #                                         base_name='parametric/', save_txt=TRUE, saveRDA=TRUE){
  #   
  #   # Create a folder to save the results
  #   dir.create(base_name, showWarnings = FALSE)
  #   
  #   
  #   # Run COMBAT
  #   # Pls run 2 options: parametric and non parametric and get the best output
  #   t1 = Sys.time()
  #   combat_output = ComBat(dat=as.matrix(myData), batch=batch, 
  #                          mod=mod, par.prior=par_prior, mean.only=mean_only, prior.plots=prior_plots)
  #   t2 = Sys.time()
  #   
  #   if(save_txt){
  #     write.table(combat_output, file = paste0(base_name,"combat_normalized_matrix.txt"), row.names = T, col.names = T, sep="\t")  
  #     print("Data saved as txt format")
  #   }
  #   if(saveRDA){
  #     save(combat_output, file = paste0(base_name,"combat_normalized.rda"))
  #     print("Data saved as rda format")
  #   }
  #   
  #   
  #   ######################
  #   ####  Execution time 
  #   ######################
  #   # Export execution time, pls load this function from evaluation_utils.R first 
  #   labeluc <- "combat_time_execution"
  #   runtime_export_func(t1, t2, labeluc, base_name)
  #   return(combat_output)
  # }
  
  
# First option, opt1 <- 'parametric/'
# mod=NULL, par.prior=TRUE, prior.plots=FALSE
combat_output <- combat_batch_effect_removal(myFilteredData, batch=mySample$batch, 
                                             mod=NULL, par_prior=TRUE, prior_plots=FALSE, mean_only=FALSE, 
                                             base_name=opt1, save_txt=TRUE, saveRDS=TRUE)

# Second option, opt2 <- 'non-parametric/'
# mod=NULL, par.prior=FALSE,prior.plots=FALSE, mean.only=TRUE
# combat_output <- combat_batch_effect_removal(pbmc@scale.data, batch=pbmc@meta.data[cells_use,'batch'],
#                             mod=NULL, par_prior=FALSE, prior_plots=FALSE, mean_only=TRUE,
#                             base_name=opt2, save_txt=FALSE, saveRDA=TRUE)

  
base_name <- opt1
# base_name <- opt2
  
  
# Create Seurat object to compute tsne, umap, pca
combat_srt <- CreateSeuratObject(raw.data = combat_output, project = "dataset4_combat")
dim(combat_srt@raw.data)  # 15558  14767
cells_use <- colnames(combat_srt@raw.data)
  
combat_srt@meta.data$batchlb <- mySample$batch
combat_srt@meta.data$celltype <- mySample$celltype
  
max(combat_srt@raw.data)
min(combat_srt@raw.data)
combat_srt@scale.data <- combat_srt@raw.data
# combat_srt <- NormalizeData(combat_srt, normalization.method = "LogNormalize")
# combat_srt <- ScaleData(object = combat_srt, vars.to.regress = c("nUMI"))
# combat_srt <- FindVariableGenes(object = combat_srt, mean.function = ExpMean, dispersion.function = LogVMR, 
                                  # x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5) #y.cutoff = 0.5
# length(combat_srt@var.genes)
  
# pc.genes = combat_srt@var.genes equal null in this scenario
  
combat_srt <- RunPCA(object = combat_srt, pcs.compute = 20, pc.genes = rownames(combat_srt@raw.data),
                     do.print = T, pcs.print = 1:5, genes.print = 5)  
  
# Visualization
combat_srt <- RunTSNE(combat_srt, reduction.use = "pca", dims.use = 1:20,perplexity=30) #do.fast = T
png(paste(base_name,"tsne_combat_parametric_v2.png",sep=""),width = 2*800, height = 2*500, res = 2*72)
p1 <- TSNEPlot(combat_srt, do.return = T, pt.size = 0.5, group.by = "batchlb")
p2 <- TSNEPlot(combat_srt, do.return = T, pt.size = 0.5, group.by = "celltype")
plot_grid(p1, p2)
dev.off()
  
combat_srt <- RunUMAP(combat_srt, reduction.use = "pca", dims.use = 1:20)
png(paste(base_name,"umap_combat_parametric_v2.png",sep=""),width = 2*800, height = 2*500, res = 2*72)
p11 <- DimPlot(object = combat_srt, reduction.use = 'umap', group.by ="batchlb")
p12 <- DimPlot(object = combat_srt, reduction.use = 'umap', group.by ="celltype")
plot_grid(p11, p12)
dev.off()
  
# Save data, easy to load and evaluate next time
save(combat_srt, file = paste0(base_name,"combat_seurat_objs.rda"))  