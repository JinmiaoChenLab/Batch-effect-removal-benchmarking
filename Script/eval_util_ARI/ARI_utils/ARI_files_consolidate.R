# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: Fourth and last function to be called in ARI pipeline
#          Consolidate all raw data and write them in a text file

#assume all ARI files are within an 'output' folder specific to a dataset

###start from within said folder
#setwd("./ARI_dataset1/ARISampled_CT/")

##eg. combinedresult<-ari_consolidate()

ari_consolidate<-function(){
  baseName<-rev(unlist(strsplit(getwd(), split = "/", fixed = TRUE)))[1]
  wds<-rev(unlist(strsplit(getwd(), split = "/", fixed = TRUE)))[2]
  filenames<-list.files(pattern = "_ARI.txt")
  bigone<-data.frame(NULL)

  for (x in 1:length(filenames)){
    temp<-read.table(file = filenames[x], header = TRUE, stringsAsFactors = FALSE)
    temp2<-temp[1:(dim(temp)[1]-1),]
    bigone<-rbind(bigone, temp2)
  }
  
  write.table(bigone, file = paste0("AllValues_", baseName, wds, ".txt"), row.names = FALSE, 
              col.names = TRUE, quote = FALSE, sep="\t")
  return(bigone)
}

