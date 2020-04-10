# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: Third function to be called in ARI pipeline
#          Following calculation of ARI scores for all batch-
#          correction methods, this function is used to 
#          produce F1 score based on normalised ARI scores 
#          Returns a CSV file containing median and F1 scores
#          for all batch-correction methods

conclude_ARISampled <- function(dir_this, plot_title){
  library(ggplot2)
  
  setwd(dir_this)
  filesneed<-list.files(pattern = "ARISampled_")
  
  method<-vector()
  medianARIbatch<-vector()
  medianARIcelltype<-vector()
  
  for (x in 1:length(filesneed)){
    temp<-read.table(filesneed[x], header = TRUE, stringsAsFactors = FALSE)
    method[x]<-as.character(temp$use_case[1])
    medianARIbatch[x]<-temp$ari_batch[21]
    medianARIcelltype[x]<-temp$ari_celltype[21]
    rm(temp)
  }
  
  # normalise values to 0 - 1
  min_batch <- min(medianARIbatch)
  max_batch <- max(medianARIbatch)
  min_cell <- min(medianARIcelltype)
  max_cell <- max(medianARIcelltype)
  medianARIbatch_norm <- (medianARIbatch-min_batch)/(max_batch-min_batch)
  medianARIcelltype_norm <- (medianARIcelltype-min_cell)/(max_cell-min_cell)
  
  # produce final fscore ARI, similar to scMerge paper
  medianfscoreARI <- (2 * (1 - medianARIbatch_norm)*(medianARIcelltype_norm))/
                          (1 - medianARIbatch_norm + medianARIcelltype_norm)
  
  sum_xy<-medianARIcelltype_norm+(1-medianARIbatch_norm)
  
  finaldf<-data.frame("ARIMethod" = method, 
                      "ARIbatchMedian" = medianARIbatch, 
                      "ARIbatchMedian_norm" = medianARIbatch_norm, 
                      "ARIcelltypeMedian" = medianARIcelltype, 
                      "ARIcelltypeMedian_norm" = medianARIcelltype_norm, 
                      "ARI_fscore" = medianfscoreARI,
                      "ARI_summedXY" = sum_xy)
  
  write.csv(finaldf, file = paste0("ARI_Sampled", plot_title, "_allmethod.csv"), row.names = FALSE)
  return(finaldf)
}
