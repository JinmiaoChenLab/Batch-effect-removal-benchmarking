# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: Second function to be called in ARI pipeline
#          Performs subsampling, extracts common cells for batch score
#          calculation, calculates the ARI scores for batch and cell type
#          Returns individual text files for each batch-correction method
#          (Specific script for Dataset 9)

####################################
# ARI Adjusted Rand Index (sampled)
####################################
# Input: myData: 20 PCs and batch, celltype columns
# cpcs: vector containing column number of PCs in myData

ari_calcul_sampled_9 <- function(myData, cpcs, isOptimal=FALSE, 
                                 method_use='resnet',
                                 base_name='', maxiter=30, 
					             celltypelb='celltype', batchlb='batch' )
{
  library(NbClust)
  library(mclust)
  set.seed(0)
    
  data9clus<-2  #best number of cluster
  
  # run function 20 times, each time extract 80% of data
  nbiters <- 20
  percent_extract <- 0.8
  
  it <- c()
  total_ari_batch <- c()

  # start loop for 20 times
  for(i in 1:nbiters) {
    
	# select cells for the subsampled dataset
    cells_extract<-sample(c(1:dim(myData)[1]), size=round(dim(myData)[1]*percent_extract), 
                          replace = F)
    cells_extract<-sort(cells_extract)
    
	# create the subsampled dataset
    myPCAExt <- myData[cells_extract,]
    
    ###############################
    # Clustering
    ###############################
    
    if(!isOptimal){  # isOptimal==FALSE
      # nbct : k equal number of unique cell types in the dataset
      clustering_result <- kmeans(x = myPCAExt[,cpcs], centers=data9clus, iter.max = maxiter)
      myPCAExt$clusterlb <- clustering_result$cluster
      
    } else if(isOptimal){
      nbclust_result<-NbClust(data=myPCAExt[,cpcs], method = "kmeans", 
	                          min.nc = max(data9clus-2, 2), max.nc = nbct+4)
      myPCAExt$clusterlb <-nbclust_result$Best.partition	  
    }
    
	# assign the current myPCAExt to a unique object so that it can be stored later 
    assign(paste0("myPCAExt",i), myPCAExt)
    print('Nb clusters: ')
    print(length(unique(myPCAExt$clusterlb)))
    
    ###############################
    # ARI
    ###############################
    
    # run ARI
    ari_batch <- mclust::adjustedRandIndex(myPCAExt[,batchlb], myPCAExt$clusterlb)
   #ari_celltype<-mclust::adjustedRandIndex(myPCAExt[,celltypelb], myPCAExt$clusterlb)
    
    it <- c(it,i)
    total_ari_batch <- c(total_ari_batch,ari_batch)
  }   # End of loop   
  
  
  # once looped 20 times to produce a total of 40 scores, next step follows:
  it <- c(it,nbiters+1)
  total_ari_batch <- c(total_ari_batch, median(total_ari_batch))

  methods <- rep(method_use, nbiters)
  methods <- c(methods,paste0(method_use,'_median'))
  
  # create final dataframe containing raw and median ARI scores
  myARI <- data.frame("use_case"=methods, 
                      "iteration"=it,
                      "ari_batch"=total_ari_batch)

  # write final dataframe to a text file
  write.table(myARI, file = paste0(base_name,method_use,"_ARI.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  
  print('Save output in folder')
  print(base_name)
  
  return(list(myARI, myPCAExt1, myPCAExt2, myPCAExt3, myPCAExt4, myPCAExt5, myPCAExt6,
              myPCAExt7, myPCAExt8, myPCAExt9, myPCAExt10, myPCAExt11, 
              myPCAExt12,  myPCAExt13, myPCAExt14, myPCAExt15, myPCAExt16,  
              myPCAExt17, myPCAExt18, myPCAExt19, myPCAExt20))
}
		 
