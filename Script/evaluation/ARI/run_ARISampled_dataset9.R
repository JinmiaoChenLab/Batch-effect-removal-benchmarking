# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: Primary script to run ARI pipeline

# clear workspace
rm(list=ls())

# source relevant functions
source("./ARI_utils/run_ARISampled_dat9.R")
source("./ARI_utils/ari_calcul_sampled_dat9.R")
source("./ARI_utils/conclude_ARISampled_dat9.R")

##############################
############ Set arguments

eval_metric <- 'ARISampled'

# read in files from this_dir
this_dir<-"./generate_PCA_tSNE_UMAP/dataset9/"

# send output to out_dir
out_dir<-"./ARI_dataset9/"

# create relevant folder within the out_dir
dir.create(paste0(out_dir, eval_metric, "_CT"),showWarnings = FALSE)
#dir.create(paste0(out_dir, eval_metric, "_OP"),showWarnings = FALSE)

##############################
############ Run 'run_ARISampled' function on relevant files 

methods_use <- 'ComBat'
fn <- 'dataset9_combat_pca.csv'
Rcombat<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'MNN_Correct'
fn <- 'dataset9_classicMNN_pca.csv'
Rmnncorrect<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'fastMNN'
fn <- 'dataset9_fastMNN_pca.csv'
Rfastmnn<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'Harmony'
fn <- 'dataset9_harmony_pca.csv'
Rharmony<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'LIGER'
fn <- 'dataset9_liger_pca.csv'
Rliger<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'limma'
fn <- 'dataset9_limma_pca.csv'
Rlimma<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'Raw'
fn <- 'dataset9_raw_pca.csv'
Rraw<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'MMD-ResNet'
fn <- 'dataset9_resnet_pca.csv'
Rresnet<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'Scanorama'
fn <- 'dataset9_scanorama_pca.csv'
Rscanorama<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'scGen'
fn <- 'dataset9_scgen_pca.csv'
Rscgen<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'scMerge'
fn <- 'dataset9_scmerge_pca.csv'
Rscmerge<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'ZINB-WaVE'
fn <- 'dataset9_zinbwave_pca.csv'
Rzinbwave<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'Seurat_2'
fn <- 'dataset9_seurat2_pca.csv'
Rseurat2<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'Seurat_3'
fn <- 'dataset9_seurat3_pca.csv'
Rseurat3<-run_ARISampled_9(fn,this_dir,out_dir,eval_metric,methods_use)
####################################################################

##############################
############ Extracting all data from all methods in dataset 

# Reads files from dir_this 
dir_this<-paste0(out_dir, eval_metric, "_CT")
#dir_this<-paste0(out_dir, eval_metric, "_OP")

# Perform the following function to produce final CSV file
wholedf<-conclude_ARISampled_9(dir_this, "Dataset_9") 

save.image(paste0("Dataset9", "_complete.RData"))

####################################################################

##############################
############ Consolidate all raw data into one file

setwd(dir_this)
rm(list=ls())

source("./ARI_utils/ARI_files_consolidate.R")
 
ari_consolidate()

### END
