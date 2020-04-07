
########## summary confusion matrix

rm(list=ls())
this.dir <- '/home/marion/Marion/Project/Hoa_batch_normalization/simulation_dataset_V3/'
setwd(this.dir)

# load libraries
library(gridExtra)

main_Fscore <- function(select){
  
  # METHODS
  vect_method <- c('raw_data','seurat3','MNN','Combat','limma','scGen','scanorama2','zinb_wave','scmerge')
  
  # HVG (all/as Seurat)
  vect_HVG <- c('all','HVG')
  
  # SIMULATIONS
  vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
  
  Fscore_list <- lapply(vect_HVG,function(HVG){
    
    base_name <- paste0('Venn_diagram/')
    
    df_all <- lapply(vect_simu,function(simu){
      
      df <- sapply(vect_method,function(method){
        
        # load all gene lists and keep upregulated
        if(method=='raw_data'){
          S3 <- read.table(paste0('Seurat_DEGs/',simu,'/','S3_batch12/',method,'_',HVG,'/degs_batch12_seurat_bimod_DEG.txt'), head=T, sep='\t')
        } else {
          S3 <- read.table(paste0('Seurat_DEGs/',simu,'/','S3_batch12/','after_',method,'_',HVG,'/degs_batch12_seurat_bimod_DEG.txt'), head=T, sep='\t')
        }
        
        geneinfo <- read.table(paste0('data/',simu,'/geneinfo.txt'), head=T)
        real_de_genes_ls <- rownames(geneinfo[(geneinfo$DEFacGroup1+geneinfo$DEFacGroup2)!=2,])
        S7 <- geneinfo[real_de_genes_ls,]
        
        if(HVG=='HVG'){
          # read HVG table 
          HVG <- read.table(paste0('data/',simu,'/counts_HVG.txt'), head=T, row.names=1)
          HVG <- colnames(HVG)
          S7 <- S7[rownames(S7) %in% HVG,] 
        } 
        
        if(select=='UP'){
          S3 <- S3[S3$avg_logFC>0,]
          S7 <- S7[S7$DEFacGroup1>S7$DEFacGroup2,]
        } else if(select=='DOWN'){
          S3 <- S3[S3$avg_logFC<0,]
          S7 <- S7[S7$DEFacGroup1<S7$DEFacGroup2,]
        } else {
          stop('select UP or DOWN DEGs')
        }
        
        GT <- S7$Gene
        norm <- S3$X
        
        # table
        TP <- sum(GT %in% norm)
        FN <- length(GT)-TP
        FP <- length(norm)-TP
        
        # recall (sensitivity)
        TPR <- round(TP/(TP + FN),3)
        
        # precision (positive predictive value)
        PPV <- round(TP/(TP + FP),3)
        
        # F-score
        Fscore <- 2*((PPV*TPR)/(PPV+TPR))
        Fscore <- round(Fscore,3)
        
        # table matrix confusion 
        data <- data.frame(TP=TP,FN=FN,FP=FP,recall=TPR,precision=PPV,Fscore=Fscore,row.names = method)
        return(data)
        
      })
      
      mef <- t(df)
      rownames(mef) <- c('Raw','Seurat 3','MNN correct','Combat','limma','scGen','Scanorama','ZINB-WaVE','scMerge')
      mef <- rbind(rep('',dim(mef)[2]),mef)
      rownames(mef)[rownames(mef)==""] <- simu
      return(list(mef,df['Fscore',]))
      
    })
    
    df_all_all <- do.call(rbind,lapply(df_all,function(l){l[[1]]}))
    
    write.table(df_all_all,paste0(base_name,'summary_confusion_matrix_',HVG,'_DEGs_',select,'.txt'),sep='\t', quote=F, row.names=T, col.names=NA)
    return(df_fscore)
    
  })
  names(Fscore_list) <- vect_HVG
  return(Fscore_list)
}


Fscore_up <- main_Fscore(select='UP') # comment write.table(df_all_all)
Fscore_down <- main_Fscore(select='DOWN') # comment write.table(df_all_all)
