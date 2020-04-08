
devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))
library(Seurat)

this.dir <- '/home/marion/Marion/Project/Hoa_batch_normalization/simulation_dataset_V3/'
setwd(this.dir)

lsdir <- list.dirs('data', recursive=FALSE) 

sapply(lsdir,function(x){
  
  # read data counts and cellinfo
  counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t')
  counts <- t(counts)
  cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t')
  cellinfo <- cellinfo[colnames(counts),]
  
  pbmc <- CreateSeuratObject(raw.data = counts, project = '',min.cells = 3,min.genes = 200)
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e2)
  
  #png(paste0(x,'/variableGenes.png'),height = 480, width = 480, res = 72)
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  #dev.off()
  length(x = pbmc@var.genes)
  
  counts_HVG <- counts[rownames(counts) %in% pbmc@var.genes,]
  write.table(t(counts_HVG), file = paste0(x,'/counts_HVG.txt'), quote=FALSE, row.names = T, col.names = T, sep="\t")

})
