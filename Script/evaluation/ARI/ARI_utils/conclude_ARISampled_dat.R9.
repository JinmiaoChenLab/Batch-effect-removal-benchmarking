# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: Third function to be called in ARI pipeline
#          Following calculation of ARI scores for all batch-
#          correction methods, this function is used to 
#          produce F1 score based on normalised ARI scores 
#          Returns a CSV file containing median and F1 scores
#          for all batch-correction methods
#          (Specific script for Dataset 9)

conclude_ARISampled_9 <- function(dir_this, plot_title){
  library(ggplot2)
  
  setwd(dir_this)
  filesneed<-list.files(pattern = "ARISampled_")
  
  method<-vector()
  medianARIbatch<-vector()
  
  for (x in 1:length(filesneed)){
    temp<-read.table(filesneed[x], header = TRUE, stringsAsFactors = FALSE)
    method[x]<-as.character(temp$use_case[1])
    medianARIbatch[x]<-temp$ari_batch[21]
    rm(temp)
  }
  
  # normalise values to 0 - 1
  min_batch <- min(medianARIbatch)
  max_batch <- max(medianARIbatch)
  medianARIbatch_norm <- (medianARIbatch-min_batch)/(max_batch-min_batch)
  
  finaldf<-data.frame("ARIMethod" = method, 
                      "ARIbatchMedian" = medianARIbatch, 
                      "ARIbatchMedian_norm" = medianARIbatch_norm) 
  
  write.csv(finaldf, file = paste0("ARI_Sampled", plot_title, "_allmethod.csv"), row.names = FALSE)
  return(finaldf)
}
