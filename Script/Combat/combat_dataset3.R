# #install.packages('devtools')
# devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))
detach("package:ezTools", unload=TRUE)
library("Seurat", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")
library(sva)

this.dir <- '/home/marion/Marion/Project/Hoa_batch_normalization/simulation_dataset_V3/'
setwd(this.dir)

lsdir <- list.dirs('data', recursive=FALSE) 

sapply(lsdir,function(x){

  x2 <- gsub('data/','',x)
  dir.create(paste0('demo_Combat/',x2), showWarnings = FALSE)
  
  selection <- c('HVG','all')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_Combat/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo
    if(s=='HVG'){
      counts <- read.table(paste0(x,'/counts_HVG.txt'), head=T, sep='\t')
    } else {
      counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t')
    }
    counts <- t(counts)
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t')
    cellinfo <- cellinfo[colnames(counts),]
    
    # Normalize each library to the median of the transcript counts across all cells 
    # Then, log transform expression values   
    print("Median normalizing counts and log-transforming")
    col_sums = apply(counts,2, sum)
    med_trans = median(col_sums)
    norm_counts = med_trans* scale(counts, center=FALSE, scale=col_sums)
    myFilteredData = log(norm_counts + 1)
    
    # pbmc <- CreateSeuratObject(raw.data = counts, project = '',min.cells = 0, min.genes = 0)
    # pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
    
    # # remove genes with variance equals to 0
    rv_genes <- which(apply(myFilteredData, 1, var)==0) # apply on normalized data
    rv_genes_names <- rownames(myFilteredData)[rv_genes]
    count_df <- myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
    
    # Run COMBAT
    t1 = Sys.time()
    combat_output = ComBat(dat=as.matrix(count_df), 
                           batch=cellinfo$Batch, 
                           mod=NULL, par.prior=TRUE, prior.plots=FALSE, mean.only=FALSE)
    t2 = Sys.time()
    
    # save the output
    save(combat_output,file=paste0('demo_Combat/',x2,'/',s,"/output.rda"))
    write.table(combat_output, file = paste0('demo_Combat/',x2,'/',s,"/output.txt"), quote=FALSE, row.names = T, col.names = T, sep="\t")
    
    # Visualization
    combat_srt <- CreateSeuratObject(raw.data = combat_output, project = "combat", min.cells = 0, min.genes = 0)
    combat_srt <- AddMetaData(object = combat_srt, metadata = subset(cellinfo,select=c('Batch','Group')))
    combat_srt <- ScaleData(object = combat_srt)
    combat_srt <- FindVariableGenes(object = combat_srt, mean.function = ExpMean, dispersion.function = LogVMR,
                                    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.1) 
    length(combat_srt@var.genes)
    
    combat_srt <- RunPCA(object = combat_srt, genes.print = 5)
    combat_srt <- RunTSNE(combat_srt, reduction.use = "pca", dims.use = 1:20, perplexity=90)
    
    png(paste0('demo_Combat/',x2,'/',s,"/tsne.png",sep=""),width = 2*800, height = 2*500, res = 2*72)
    p1 <- TSNEPlot(combat_srt, do.return = T, pt.size = 0.5, group.by = "Batch")
    p2 <- TSNEPlot(combat_srt, do.return = T, pt.size = 0.5, group.by = "Group")
    print(plot_grid(p1, p2))
    dev.off()
    
    tsne_df <- data.frame(combat_srt@dr$tsne@cell.embeddings)
    tsne_df$batch <- combat_srt@meta.data$Batch
    tsne_df$celltype <- combat_srt@meta.data$Group
    write.table(tsne_df, file = paste0('demo_Combat/',x2,'/',s,"/tsne.txt"), row.names = T, col.names = T, sep="\t")
    
  })
  
})
