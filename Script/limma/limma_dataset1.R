## Title : limma dataset 1
## Author : Michelle Goh
## Date : 28/08/2019


### Clear environment
rm(list=ls())


### Packages
require(edgeR)
require(Seurat)
packageVersion('Seurat')


### Setup directories
this.dir <- '/home/hoa/hoatran/demo_normalization/michelle/demo_limma/'
setwd(this.dir)
data_dir <- '/home/hoa/hoatran/demo_normalization/dataset/dataset1_uc3/'
base_name <- 'limma_dataset1/'
dir.create(base_name, showWarnings = FALSE)


### load limma functions files
utils_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/michelle/'
source(paste0(utils_dir,'limma_functions.R'))


### Load files
TPM_file <- 'dataset1_sm_uc3.txt'  # replace by link to dataset
sample_file <- 'sample_sm_uc3.txt' # replace by link to dataset
myData <- read.table(paste0(data_dir,TPM_file),sep="\t",header=T,row.names=1,check.names=F)
mySample <- read.table(paste0(data_dir,sample_file),sep="\t",header=T,row.names=1,check.names=F)
# mySample$celltype <- mySample$ct


### Filter, norm & log data then return a DGElist object
dge <- DGEList(counts=myData)
group <- mySample$batch
design <- model.matrix(~group) 
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method = "TMMwzp")
v <- voom(dge, design, plot=TRUE)


### limma function 
limma <- limma_batch_effect_removal(v, mySample, base_name, save_txt=FALSE, saveRDS=TRUE)


### Create Seurat object to compute tsne, umap, pca
limma_srt <- CreateSeuratObject(raw.data = limma, project = "dataset1_limma")
dim(limma_srt@raw.data)  # 2562 576
cells_use <- colnames(limma_srt@raw.data)

limma_srt@meta.data$batchlb <- mySample$batch
limma_srt@meta.data$celltype <- mySample$celltype 

max(limma_srt@raw.data)
min(limma_srt@raw.data)
limma_srt@scale.data <- limma_srt@raw.data

limma_srt <- RunPCA(object = limma_srt, pcs.compute = 20, pc.genes = rownames(limma_srt@raw.data),
                     do.print = T, pcs.print = 1:5, genes.print = 5)  

# Visualization
limma_srt <- RunTSNE(limma_srt, reduction.use = "pca", dims.use = 1:20,perplexity=30) #do.fast = T
png(paste(base_name,"tsne_limma_voom.png",sep=""),width = 2*800, height = 2*500, res = 2*72)
p1 <- TSNEPlot(limma_srt, do.return = T, pt.size = 0.5, group.by = "batchlb")
p2 <- TSNEPlot(limma_srt, do.return = T, pt.size = 0.5, group.by = "celltype")
plot_grid(p1, p2)
dev.off()

limma_srt <- RunUMAP(limma_srt, reduction.use = "pca", dims.use = 1:20)
png(paste(base_name,"umap_limma_voom.png",sep=""),width = 2*800, height = 2*500, res = 2*72)
p11 <- DimPlot(object = limma_srt, reduction.use = 'umap', group.by ="batchlb")
p12 <- DimPlot(object = limma_srt, reduction.use = 'umap', group.by ="celltype")
plot_grid(p11, p12)
dev.off()

# Save data, easy to load and evaluate next time
save(limma_srt, file = paste0(base_name,"limma_seurat_objs.rda"))  