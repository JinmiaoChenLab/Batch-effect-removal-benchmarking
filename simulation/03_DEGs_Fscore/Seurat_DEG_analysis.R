
seurat_analysis_deg <- function(TPM=TPM,
                            sample=sample,
                            group_col="Site",
                            base_name="GC_cluster7", 
                            mt="^MT-", 
                            scale=F, 
                            scale_factor=1e4,
                            group1=NULL,
                            group2=NULL,
                            test.use="bimod",
                            col.low = "#FF00FF",
                            col.mid = "#000000", col.high = "#FFFF00",
                            logfc.threshold = 0.1,
                            plotDEGhm=T,topn=NULL,
                            raw_data=F){

  library(Seurat)
  library(dplyr)
  library(Matrix)

  pbmc.data <- TPM
  
  common_cell <- intersect(row.names(sample),colnames(pbmc.data))
  pbmc.data <- pbmc.data[,common_cell]
  if(ncol(sample)==1) sample <- cbind(sample,sample)
  sample <- sample[common_cell,]

  if((!is.null(group1)) & (!is.null(group2))){
    group12 <- row.names(sample)[as.character(sample[,group_col])==as.character(group1) | as.character(sample[,group_col])==as.character(group2)]
    pbmc.data <- pbmc.data[,group12]
  }
  
  if(!is.null(sample)){
    colnames(pbmc.data) <- paste(sample[colnames(pbmc.data),group_col],colnames(pbmc.data),sep="_")
  }
  
  pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 0, min.genes = 0, project = base_name)
  
  if(raw_data=T){
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = scale_factor)
  } else {
    # do no normalize or scale the output from batch effect correction methods (already done inside the methods)
    pbmc@data <- pbmc@raw.data
  }
  
  if(is.null(group1)&is.null(group2)){
     pbmc.markers <- FindAllMarkers(object = pbmc, test.use=test.use, logfc.threshold = logfc.threshold)
  }else{
	   pbmc.markers <- FindMarkers(object = pbmc, ident.1 = group1, ident.2 = group2,test.use=test.use,logfc.threshold = logfc.threshold)
	}

  pbmc.markers <- pbmc.markers[!pbmc.markers$avg_logFC=='Inf',]
  
  if(sum(pbmc.markers$p_val_adj<0.05,na.rm = T)==0){
     print("No DEGs with adjusted p values < 0.05")
     return(0)
  }
  ####################

  #rotate = FALSE

  if(is.null(group1) & is.null(group2)){

    if(length(levels(pbmc.markers$cluster))>2){
      deg <- pbmc.markers[pbmc.markers$p_val_adj<0.05 & pbmc.markers$avg_logFC>0,]
      rotate=TRUE
    } else {
      deg <- pbmc.markers[pbmc.markers$p_val_adj<0.05,]
    }
    deg <- deg[order(deg$cluster,deg$avg_logFC,decreasing = T),]
    genelist <- unique(deg$gene)

    # keep topn
    if(!is.null(topn)){
      degtop <- deg[order(deg$avg_logFC,decreasing = T),]
      degtop <- degtop[!duplicated(degtop$gene),]
      topnn <- topn
      if(topn>dim(degtop)[1]){topnn <- dim(degtop)[1]}
      degtop <- degtop[1:topnn,]
      degtop <- degtop[order(degtop$cluster,degtop$avg_logFC,decreasing = T),]
      top50deg <- unique(degtop$gene)
    }

  } else {

    deg <- pbmc.markers[pbmc.markers$p_val_adj<0.05,]
    deg <- deg[order(deg$avg_logFC,decreasing = T),]
    genelist <- rownames(deg)

    if(!is.null(topn)){
      degtop <- deg[order(abs(deg$avg_logFC),decreasing = T),]
      topnn <- topn
      if(topn>dim(degtop)[1]){topnn <- dim(degtop)[1]}
      degtop <- degtop[1:topnn,]
      degtop <- degtop[order(degtop$avg_logFC,decreasing = T),]
      top50deg <- rownames(degtop)

    }
  }
  
  write.table(deg,paste(base_name,"_seurat_bimod_DEG.txt",sep=""),sep="\t",col.names=NA, quote=F)
  
}
