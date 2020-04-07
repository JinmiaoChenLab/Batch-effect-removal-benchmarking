
rm(list=ls())

library(SingleCellExperiment)
library(scater)
#BiocManager::install("scMerge")
library(scMerge)

this.dir <- '/scbio4/home/marion/Hoa_batch_normalization/simulation_dataset_V3/'
setwd(this.dir)

lsdir <- list.dirs('data', recursive=FALSE) 

sapply(lsdir,function(x){
  
  x2 <- gsub('data/','',x)
  dir.create(paste0('demo_scmerge/',x2), showWarnings = FALSE)
  
  selection <- c('HVG','all')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_scmerge/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo
    if(s=='HVG'){
      counts <- read.table(paste0(x,'/counts_HVG.txt'), head=T, sep='\t')
    } else {
      counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t')
      counts_HVG <- read.table(paste0(x,'/counts_HVG.txt'), head=T, sep='\t')
    }
    counts <- t(counts)
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t')
    cellinfo <- cellinfo[colnames(counts),]
    cellinfo <- subset(cellinfo,select=c('Batch','Group'))
    names(cellinfo) <- c('batch','cell_type')
    
    counts <- counts[rowSums(counts)>1,]
    
    # remove genes with variance equals to 0
    rv_genes <- which(apply(counts, 1, var)==0) 
    rv_genes_names <- rownames(counts)[rv_genes]
    count_df <- counts[!(rownames(counts) %in% rv_genes_names),]
    
    count_df <- count_df[,rownames(cellinfo)]
    
    sce <- SingleCellExperiment(assays = list(counts = count_df),colData = cellinfo)
    sce <- normalize(sce)
    
    counts_b1 <- count_df[,colnames(count_df) %in% rownames(cellinfo[cellinfo$batch=='Batch1',])]
    counts_b2 <- count_df[,colnames(count_df) %in% rownames(cellinfo[cellinfo$batch=='Batch2',])]
    data <- list(counts_b1,counts_b2)
    
    geneinfo <- read.table(paste0(x,'/geneinfo.txt'), head=T, sep='\t')
    geneinfo_wobatch <- geneinfo[geneinfo$DEFacGroup1==1 & geneinfo$DEFacGroup2==1,]
    counts2 <- count_df[rownames(count_df) %in% as.character(geneinfo_wobatch$Gene),]
    gene_lowvar <- names(sort(apply(counts2,1,var), decreasing=F)[1:200])
    
    if(s=='HVG'){
      # Run ScMerge
      t1 = Sys.time()
      scMerge_res <- scMerge(
        sce_combine = sce, 
        ctl = gene_lowvar,
        kmeansK = c(2,2),
        assay_name = "scMerge_res",
        replicate_prop = 0.5, verbose=T,
        marker = rownames(counts))
      t2 = Sys.time()
      print(t2-t1)
    } else {
      # Run ScMerge
      t1 = Sys.time()
      scMerge_res <- scMerge(
        sce_combine = sce, 
        ctl = gene_lowvar,
        kmeansK = c(2,2),
        assay_name = "scMerge_res",
        replicate_prop = 0.5, verbose=T,
        marker = colnames(counts_HVG))
      t2 = Sys.time()
      print(t2-t1)
    }
    
    # save the output
    save(scMerge_res, file = paste0('demo_scmerge/',x2,'/',s,"/output.rda"))
    scmerge_norm <- scMerge_res@assays[["scMerge_res"]]
    write.table(scmerge_norm, file = paste0('demo_scmerge/',x2,'/',s,"/output.txt"), row.names = T, col.names = T, sep="\t")
    
    # Visualization
    scMerge_res = scater::runTSNE(scMerge_res, exprs_values = "scMerge_res")
    df <- data.frame(scMerge_res@reducedDims$TSNE)
    colnames(df) <- paste0('tSNE_',c(1:2))
    rownames(df) <- rownames(scMerge_res@colData)
    df <- ezTools::ezcbind(df,scMerge_res@colData)
    write.table(df, file = paste0('demo_scmerge/',x2,'/',s,"/tsne.txt"), row.names = T, col.names = T, sep="\t")
    
    png(paste0('demo_scmerge/',x2,'/',s,"/tsne.png",sep=""),width = 2*800, height = 2*500, res = 2*72, type='cairo')
    p1 <- ggplot(df, aes(x=tSNE_1,y=tSNE_2, color=batch)) + geom_point()
    p2 <- ggplot(df, aes(x=tSNE_1,y=tSNE_2, color=cell_type)) + geom_point()
    print(plot_grid(p1, p2))
    dev.off()
    
  })
  
})


