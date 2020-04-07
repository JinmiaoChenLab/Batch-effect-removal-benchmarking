
rm(list=ls())

library(limma)
# devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))
library(Seurat)
packageVersion('Seurat')

this.dir <- '/acrc/jinmiao/CJM_lab/Marion/Project/Hoa_batch_normalization/simulation_dataset_V3/'
setwd(this.dir)

lsdir <- list.dirs('data', recursive=FALSE) 

sapply(lsdir,function(x){
  
  x2 <- gsub('data/','',x)
  dir.create(paste0('demo_limma/',x2), showWarnings = FALSE)
  
  selection <- c('HVG','all')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_limma/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo
    if(s=='HVG'){
      counts <- read.table(paste0(x,'/counts_HVG.txt'), head=T, sep='\t')
    } else {
      counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t')
    }
    counts <- t(counts)
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t')
    cellinfo <- cellinfo[colnames(counts),]
    cellinfo <- subset(cellinfo,select=c('Batch','Group'))
    names(cellinfo) <- c('batch','cell_type')
    
    # input limma : normalized matrix
    pbmc <- CreateSeuratObject(raw.data = counts, meta.data = cellinfo, project = "Limma", min.cells = 0, min.genes = 0)
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
    pbmc <- ScaleData(object = pbmc)
    output_df <- pbmc@data
    
    # run limma
    lm_df <- limma::removeBatchEffect(as.matrix(output_df), factor(cellinfo$batch))

    # save the output
    limma_srt <- CreateSeuratObject(raw.data = lm_df, meta.data = cellinfo, project = "Limma_normalized")
    save(limma_srt, file = paste0('demo_limma/',x2,'/',s,"/output.rda"))
    write.table(lm_df, file = paste0('demo_limma/',x2,'/',s,"/output.txt"), row.names = T, col.names = T, sep="\t")
    
    # Visualization
    limma_srt@data <- limma_srt@raw.data
    limma_srt@scale.data <- limma_srt@raw.data
    limma_srt <- FindVariableGenes(object = limma_srt, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.1)
    length(limma_srt@var.genes)
    limma_srt <- RunPCA(object = limma_srt, pcs.compute = 20, pc.genes = limma_srt@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5) 
    limma_srt <- RunTSNE(limma_srt, reduction.use = "pca", dims.use = 1:20, do.fast = T, perplexity=30)
    
    png(paste0('demo_limma/',x2,'/',s,"/tsne.png",sep=""),width = 2*800, height = 2*500, res = 2*72)
    p1 <- TSNEPlot(limma_srt, do.return = T, pt.size = 0.5, group.by = "batch")
    p2 <- TSNEPlot(limma_srt, do.return = T, pt.size = 0.5, group.by = "cell_type")
    print(plot_grid(p1, p2))
    dev.off()
    
    tsne_df <- data.frame(limma_srt@dr$tsne@cell.embeddings)
    tsne_df$batch <- limma_srt@meta.data$batch
    tsne_df$celltype <- limma_srt@meta.data$cell_type
    write.table(tsne_df, file = paste0('demo_limma/',x2,'/',s,"/tsne.txt"), row.names = T, col.names = T, sep="\t")
  })
  
})