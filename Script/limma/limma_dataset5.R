  # Limma batch effect removal 
  # Input: log normalize matrix 
  # Documentation: https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/removeBatchEffect
  # Hoa Tran
  
  
  library(limma)
  library(Seurat)
  packageVersion('Seurat')
  
  # Set working directory
  rm(list=ls())
  this.dir <- '/home/hoa/hoatran/demo_normalization/michelle/'
  setwd(this.dir)
  
  source('limma_functions.R')  
  
  # create base directory to save the results 
  base_name <- 'results_dataset5_human_pbmc/' 
  dir.create(base_name, showWarnings = FALSE)
  


  data_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/dataset/dataset5_human_pbmc/final/'
  
  sample_b1 <- read.table(paste0(data_dir,'b1_celltype.txt'), sep = "\t", row.names = 1, head=T, check.names = FALSE) 
  # dim(sample_b1)
  sample_b1$batch <- 'batch1'
  sample_b2 <- read.table(paste0(data_dir,'b2_celltype.txt'), sep = "\t", row.names = 1, head=T, check.names = FALSE)
  # dim(sample_b2)
  sample_b2$batch <- 'batch2'
  
  data_b1 <- read.table(paste0(data_dir,"b1_exprs.txt"), sep = "\t", row.names = 1, head=T, check.names = FALSE)
  # dim(data_b1)
  data_b2 <- read.table(paste0(data_dir,"b2_exprs.txt"), sep = "\t", row.names = 1, head=T, check.names = FALSE)
  # dim(data_b2)
  
  # Verification
  # sum(rownames(data_b1)==rownames(data_b2))
  # sum(colnames(sample_b1)==colnames(sample_b2))
  myData <- cbind(data_b1, data_b2)
  mySample <- rbind(sample_b1, sample_b2)
  mySample$celltype <- mySample$CellType
  mySample$batchlb <- mySample$batch
  
  
  # Filter data, return a Seurat object
  min_cells = 10
  # min_genes = 300
  min_genes = 0
  myFilteredData <- limma_filter_data(myData, mySample, min_cells, min_genes)
  mySample <- mySample[colnames(myFilteredData),]
  # dim(mySample)
  # View(head(mySample))
  limma <- limma_batch_effect_removal(myFilteredData, mySample, base_name, save_txt=FALSE, saveRDS=TRUE)
  
  # Create Seurat object to compute tsne, umap, pca
  limma_srt <- CreateSeuratObject(raw.data = limma, project = "dataset5_limma")
  dim(limma_srt@raw.data)  # 5000  15476
  cells_use <- colnames(limma_srt@raw.data)
  
  limma_srt@meta.data$batchlb <- mySample$batch
  limma_srt@meta.data$celltype <- mySample$celltype
  
  max(limma_srt@raw.data)
  min(limma_srt@raw.data)
  limma_srt@scale.data <- limma_srt@raw.data
  # combat_srt <- NormalizeData(combat_srt, normalization.method = "LogNormalize")
  # combat_srt <- ScaleData(object = combat_srt, vars.to.regress = c("nUMI"))
  # combat_srt <- FindVariableGenes(object = combat_srt, mean.function = ExpMean, dispersion.function = LogVMR, 
  # x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.1) #y.cutoff = 0.5
  # length(combat_srt@var.genes)
  
  # pc.genes = combat_srt@var.genes equal null in this scenario
  
  limma_srt <- RunPCA(object = limma_srt, pcs.compute = 20, pc.genes = rownames(limma_srt@raw.data),
                       do.print = T, pcs.print = 1:5, genes.print = 5)  
  
  # Visualization
  limma_srt <- RunTSNE(limma_srt, reduction.use = "pca", dims.use = 1:20,perplexity=30) #do.fast = T
  png(paste(base_name,"tsne_limma.png",sep=""),width = 2*800, height = 2*500, res = 2*72)
  p1 <- TSNEPlot(limma_srt, do.return = T, pt.size = 0.5, group.by = "batchlb")
  p2 <- TSNEPlot(limma_srt, do.return = T, pt.size = 0.5, group.by = "celltype")
  plot_grid(p1, p2)
  dev.off()
  
  limma_srt <- RunUMAP(limma_srt, reduction.use = "pca", dims.use = 1:20)
  png(paste(base_name,"umap_limma.png",sep=""),width = 2*800, height = 2*500, res = 2*72)
  p11 <- DimPlot(object = limma_srt, reduction.use = 'umap', group.by ="batchlb")
  p12 <- DimPlot(object = limma_srt, reduction.use = 'umap', group.by ="celltype")
  plot_grid(p11, p12)
  dev.off()
  
  # Save data, easy to load and evaluate next time
  save(limma_srt, file = paste0(base_name,"limma_seurat_objs.rda"))  