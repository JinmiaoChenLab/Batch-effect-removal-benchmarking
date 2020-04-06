rm(list=ls())

library(Seurat)
library(dplyr)

expression_matrix <- Read10X(data.dir = "C:/Users/jinmiao/Desktop/temp")
write.table(expression_matrix,"umi.txt",sep="\t")
seurat_object = CreateSeuratObject(raw.data = expression_matrix)