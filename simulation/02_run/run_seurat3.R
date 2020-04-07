
rm(list=ls())

library(Seurat)  # Seurat v3 version
library(cowplot)
library(ggplot2)

this.dir <- '/home/marion/Marion/Project/Hoa_batch_normalization/simulation_dataset_V3/'
setwd(this.dir)

lsdir <- list.dirs('data', recursive=FALSE) 

sapply(lsdir,function(x){
  
  x2 <- gsub('data/','',x)
  dir.create(paste0('demo_seurat3/',x2), showWarnings = FALSE)
  
  selection <- c('HVG','all')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_seurat3/',x2,'/',s), showWarnings = FALSE)
    
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
    
    pbmc.list <- SplitObject(pbmc, split.by = "Batch")
    pbmc.list$Batch1 <- NormalizeData(pbmc.list$Batch1, normalization.method = "LogNormalize", verbose=FALSE)
    pbmc.list$Batch1@assays$RNA@var.features <- rownames(counts)
    #pbmc.list$Batch1 <- FindVariableFeatures(pbmc.list$Batch1, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
    pbmc.list$Batch2 <- NormalizeData(pbmc.list$Batch2, normalization.method = "LogNormalize", verbose=FALSE)
    pbmc.list$Batch2@assays$RNA@var.features <- rownames(counts)
    #pbmc.list$Batch2 <- FindVariableFeatures(pbmc.list$Batch2, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
    
    # Run Seurat V3 integration
    t1 = Sys.time()
    immune.anchors <- FindIntegrationAnchors(object.list = pbmc.list, dims = 1:20)
    immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
    t2 = Sys.time()
    print(t2-t1)
    
    # save the output
    save(immune.combined, file = paste0('demo_seurat3/',x2,'/',s,"/output.rda"))
    seuratv3_integrated <- immune.combined@assays$integrated@data
    write.table(seuratv3_integrated, file = paste0('demo_seurat3/',x2,'/',s,"/output.txt"), row.names = T, col.names = T, sep="\t")
    
    # Visualization
    DefaultAssay(immune.combined) <- "integrated"
    immune.combined <- ScaleData(immune.combined, verbose = FALSE)
    immune.combined <- RunPCA(immune.combined, npcs = 20, verbose = FALSE)
    immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:20)
    
    png(paste0('demo_seurat3/',x2,'/',s,"/tsne.png",sep=""),width = 2*800, height = 2*500, res = 2*72)
    p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "Batch", pt.size = 0.5)
    p2 <- DimPlot(immune.combined, reduction = "tsne", group.by = "Group", pt.size = 0.5)
    print(plot_grid(p1, p2))
    dev.off()
    
    tsne_df <- data.frame(immune.combined@reductions$tsne@cell.embeddings)
    tsne_df$batch <- immune.combined@meta.data$Batch
    tsne_df$celltype <- immune.combined@meta.data$Group
    write.table(tsne_df, file = paste0('demo_seurat3/',x2,'/',s,"/tsne.txt"), row.names = T, col.names = T, sep="\t")
  })
  
})


