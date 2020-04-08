
rm(list=ls())

library(scran)
library(scales)
require(Rtsne)
library(Seurat)
library(ggplot2)
library(cowplot)

# BiocManager::install(version='devel')
# BiocManager::install("batchelor")
library(batchelor)

lsdir <- list.dirs('Data/dataset3/', recursive=FALSE) 

sapply(lsdir,function(x){
  
  x2 <- gsub('Data/dataset3/','',x)
  dir.create(paste0('demo_MNN/',x2), showWarnings = FALSE)

  selection <- c('HVG','all')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_MNN/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo
    if(s=='HVG'){
      counts <- read.table(paste0(x,'/counts_HVG.txt'), head=T, sep='\t')
    } else {
      counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t')
    }
    counts <- t(counts)
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t')
    cellinfo <- cellinfo[colnames(counts),]
    
    pbmc <- CreateSeuratObject(counts = counts, project = '', min.cells = 0, min.features = 0, meta.data = subset(cellinfo,select=c('Batch','Group')))
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")
    #pbmc <- ScaleData(object = pbmc)
    #pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    
    pbmc.list <- SplitObject(pbmc, split.by = "Batch")
    myData1 <- pbmc.list[[1]]@assays$RNA@data
    myData2 <- pbmc.list[[2]]@assays$RNA@data
    
    # Run MNN
    t1 = Sys.time()
    out.mnn.total <- batchelor::mnnCorrect(myData1, myData2, k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE)
    t2 = Sys.time()
    print(t2-t1)
    
    # save the output
    save(out.mnn.total, file = paste0('demo_MNN/',x2,'/',s,"/output.rda"))
    corre.mnn <- out.mnn.total@assays[['corrected']]
    write.table(corre.mnn, file = paste0('demo_MNN/',x2,'/',s,"/output.txt"), row.names = T, col.names = T, sep="\t")
    
    # Visualization
    set.seed(0)
    all.dists2.c <- as.matrix(dist(t(corre.mnn)))
    tsne.c <- Rtsne(all.dists2.c, is_distance=TRUE, perplexity = 30)  # as suggestion from MNN paper
    tsne_df <- data.frame(tSNE_1=tsne.c$Y[,1],tSNE_2=tsne.c$Y[,2])
    rownames(tsne_df) <- colnames(corre.mnn)
    tsne_df$batch <- cellinfo[rownames(tsne_df),'Batch']
    tsne_df$celltype <- cellinfo[rownames(tsne_df),'Group']
    write.table(tsne_df, file = paste0('demo_MNN/',x2,'/',s,"/tsne.txt"), row.names = T, col.names = T, sep="\t")
    
    png(paste0('demo_MNN/',x2,'/',s,"/tsne.png",sep=""),width = 2*800, height = 2*500, res = 2*72, type='cairo')
    p1 <- ggplot(tsne_df, aes(x=tSNE_1,y=tSNE_2, color=batch)) + geom_point()
    p2 <- ggplot(tsne_df, aes(x=tSNE_1,y=tSNE_2, color=celltype)) + geom_point()
    print(plot_grid(p1, p2))
    dev.off()
  })
  
})


