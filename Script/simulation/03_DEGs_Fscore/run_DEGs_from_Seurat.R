
rm(list=ls())

source('Seurat_DEG_analysis.R')

######## S3 batch12 after normalization 

# METHODS
vect_method <- c('Combat','limma','MNN','scanorama','scGen','scmerge','seurat3','zinb_wave')

# HVG (all/as Seurat)
#vect_HVG <- c('all','HVG')
vect_HVG <- c('all')

sapply(vect_method,function(method){
  sapply(vect_HVG,function(HVG){
    
    base_name <- paste0('Seurat_DEGs/')
    dir.create(base_name, showWarnings = FALSE)
    base_name <- paste0(base_name,'S3_batch12/')
    dir.create(base_name, showWarnings = FALSE)
    base_name <- paste0(base_name,'after_',method,'_',HVG,'/')
    dir.create(base_name, showWarnings = FALSE)
    
    # import data, sample
    inside <- list.files(paste0('demo_',method,'/',HVG), recursive=FALSE)
    if(grepl('.txt',inside[grep('output\\.((txt)|(csv))',inside)])){
      output_batch12 <- read.table(file = paste0('demo_',method,'/',HVG,'/output.txt'),sep="\t",header=T,row.names=1,check.names = F)
    } else {
      output_batch12 <- read.csv(file = paste0('demo_',method,'/',HVG,'/output.csv'),sep=",",header=T,row.names=1,check.names = F)
      output_batch12 <- t(output_batch12)
    }
    
    sample_batch12 <- read.table(paste0("data/cellinfo.txt"),sep="\t",header=T,row.names=1)
    
    #### DEGs Group1 vs Group2 (batch 1 + batch 2)
    seurat_analysis_deg(TPM=output_batch12,
                        sample=sample_batch12,
                        group_col="Group",
                        base_name=paste0(base_name,'degs_batch12'),  
                        group1='Group1',
                        group2='Group2',
                        test.use="bimod",
                        logfc.threshold = 0)
  })
})


######## S3 batch1+batch2 raw data

# HVG (all/as Seurat)
#vect_HVG <- c('all','HVG')
vect_HVG <- c('all')

sapply(vect_HVG,function(HVG){
  
  base_name <- paste0('Seurat_DEGs/')
  dir.create(base_name, showWarnings = FALSE)
  base_name <- paste0(base_name,'S3_batch12/')
  dir.create(base_name, showWarnings = FALSE)
  base_name <- paste0(base_name,'raw_data_',HVG,'/')
  dir.create(base_name, showWarnings = FALSE)
  
  sample <- read.table(paste0("data/cellinfo.txt"),sep="\t",header=T,row.names=1)
  
  # import data, sample
  if(HVG=='HVG'){
    data <- read.table(file = paste0('data/counts.txt'),sep="\t",header=T,row.names=1,check.names = F)
  } else {
    data <- read.table(file = paste0('data/counts_HVG.txt'),sep="\t",header=T,row.names=1,check.names = F)
  }
  data <- t(data)
  
  #### DEGs Group1 vs Group2 (batches 1+2)
  seurat_analysis_deg(TPM=data,
                      sample=sample,
                      group_col="Group",
                      base_name=paste0(base_name,'degs_batch12_'),
                      group1='Group1',
                      group2='Group2',
                      test.use="bimod",
                      logfc.threshold = 0,
                      raw_data=T)
})
