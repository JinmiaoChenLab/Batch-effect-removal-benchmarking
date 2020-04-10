# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: First function to be called in ARI pipeline
#          Reads in relevant dataset and grabs essential columns
#          Calls the following function 'ari_calcul_sampled'
#          (Specific script for Dataset 9)

run_ARISampled_9 <- function(fn, this_dir, out_dir, eval_metric, methods_use){
  setwd(this_dir)
  dataset_no<-rev(strsplit(this_dir, split = "/")[[1]])[1]
  
  thisData <- read.csv(paste0(this_dir,'/',fn), head=T, row.names = 1, check.names = FALSE)
  
  # Get relevant columns
  colPCA <- grep('([Pp][Cc]_?)|(V)|(Harmony)|(W)',colnames(thisData))
  colPCA <- colPCA[1:20]
  
  colnames(thisData)[grep('[cC]ell_?[tT]ype',colnames(thisData))] <- 'celltype'
  colnames(thisData)[grep('([bB]atch)|(BATCH)|(batchlb)',colnames(thisData))] <- 'batch'
  
  #setwd(paste0(out_dir, '/', eval_metric, "_OP"))
  #temp<-ari_calcul_sampled_9(myData=thisData, cpcs=colPCA, isOptimal=TRUE,
  #                           method_use = methods_use,  
  #                           base_name=paste0(dataset_no, eval_metric, '_OP_'))
  setwd(paste0(out_dir, '/', eval_metric, "_CT"))
  temp<-ari_calcul_sampled_9(myData=thisData, cpcs=colPCA, isOptimal=FALSE, 
                             method_use = methods_use,  
                             base_name=paste0(dataset_no, eval_metric, '_CT_'))
  return(temp)
}
