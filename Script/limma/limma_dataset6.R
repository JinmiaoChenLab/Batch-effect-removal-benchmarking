# Batch merge #
setwd("/home/hoa/hoatran/demo_normalization/michelle/demo_combat/")

read_dir = "/home/hoa/hoatran/demo_normalization/dataset/dataset12_cell_line/final/"
b1_exprs_file = "b1_exprs.txt"
b2_exprs_file = "b2_exprs.txt"
b3_exprs_file = "b3_exprs.txt"
b1_celltype_file = "b1_celltype.txt"
b2_celltype_file = "b2_celltype.txt"
b3_celltype_file = "b3_celltype.txt"

b1_exprs <- read.table(paste0(read_dir,b1_exprs_file), sep="\t", header=T, row.names=1, check.names=F)
b2_exprs <- read.table(paste0(read_dir,b2_exprs_file), sep="\t", header=T, row.names=1, check.names=F)
b3_exprs <- read.table(paste0(read_dir,b3_exprs_file), sep="\t", header=T, row.names=1, check.names=F)
b1_celltype <- read.table(paste0(read_dir,b1_celltype_file), sep="\t", header=T, row.names=1, check.names=F)
b2_celltype <- read.table(paste0(read_dir,b2_celltype_file), sep="\t", header=T, row.names=1, check.names=F)
b3_celltype <- read.table(paste0(read_dir,b3_celltype_file), sep="\t", header=T, row.names=1, check.names=F)

b1_batch_label <- rep('Batch_1', ncol(b1_exprs))
b1_batch_label <- data.frame("batch"=b1_batch_label)
row.names(b1_batch_label) <- colnames(b1_exprs)
b1_metadata = cbind(b1_celltype, b1_batch_label)

b2_batch_label <- rep('Batch_2', ncol(b2_exprs))
b2_batch_label <- data.frame("batch"=b2_batch_label)
row.names(b2_batch_label) <- colnames(b2_exprs)
b2_metadata = cbind(b2_celltype, b2_batch_label)

b3_batch_label <- rep('Batch_3', ncol(b3_exprs))
b3_batch_label <- data.frame("batch"=b3_batch_label)
row.names(b3_batch_label) <- colnames(b3_exprs)
b3_metadata = cbind(b3_celltype, b3_batch_label)

# expr_mat = cbind(b1_exprs,b2_exprs,b3_exprs)
#batchlabel_v = c(b1_batch_label,b2_batch_label,b3_batch_label)
#celltype_v = c(b1_celltype,b2_celltype,b3_celltype)
# data_meta = rbind(b1_metadata, b2_metadata, b3_metadata)


## Title : limma dataset 6
## Author : Michelle Goh
## Date : 28/08/2019


### Packages
require(edgeR)
require(Seurat)
packageVersion('Seurat')


### Setup directories
this.dir <- '/home/hoa/hoatran/demo_normalization/michelle/demo_limma/'
setwd(this.dir)
data_dir <- '/home/hoa/hoatran/demo_normalization/dataset/dataset12_cell_line/final/'
base_name <- 'limma_dataset6/'
dir.create(base_name, showWarnings = FALSE)


### load limma functions files
utils_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/michelle/'
source(paste0(utils_dir,'limma_functions.R'))


### Load files
# TPM_file <- 'myData_pancreatic_5batches.txt'  # replace by link to dataset
# sample_file <- 'mySample_pancreatic_5batches.txt' # replace by link to dataset
myData <- cbind(b1_exprs,b2_exprs,b3_exprs)
mySample <- rbind(b1_metadata, b2_metadata, b3_metadata)
mySample$celltype <- mySample$CellType


### Filter, norm & log data then return a DGElist object
dge <- DGEList(counts=myData)
group <- mySample$batch
design <- model.matrix(~group)
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method = "TMMwzp")
v <- voom(dge, design, plot=TRUE)


### limma function 
limma <- limma_batch_effect_removal(v, mySample, base_name, save_txt=FALSE, saveRDS=TRUE)


### Create Seurat object to compute tsne, umap, pca
limma_srt <- CreateSeuratObject(raw.data = limma, project = "dataset6_limma")
dim(limma_srt@raw.data)  # 146 9531
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