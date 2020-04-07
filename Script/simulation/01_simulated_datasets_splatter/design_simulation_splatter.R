
# Tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splat_params.R

rm(list=ls())

library(splatter)

this.dir <- "/acrc/jinmiao/CJM_lab/Marion/Project/Hoa_batch_normalization/simulation_dataset_V3/data/"
setwd(this.dir)

# Get all datasets in dataset3_simulation_v2/ folders
# Ex: dataset3_simulation_v2/ 
# Ex: dataset3_simulation_v2/simul1_dropout_005_b1_500_b2_900/    unbalanced number of cells in 2 batches, small dropout
# Ex: dataset3_simulation_v2/simul2_dropout_025_b1_500_b2_900/    unbalanced number of cells in 2 batches, large dropout
# Ex: dataset3_simulation_v2/simul3_dropout_005_b1_500_b2_450/    balanced number of cells in 2 batches, small dropout
# Ex: dataset3_simulation_v2/simul3_dropout_025_b1_500_b2_450/    unbalanced number of cells in 2 batches, large dropout
# Ex: dataset3_simulation_v2/simul4_dropout_005_b1_80_b2_400/     small number of cells in 2 batches, small dropout
# Ex: dataset3_simulation_v2/simul4_dropout_025_b1_80_b2_400/     small number of cells in 2 batches, large dropout

simulate <- function(nGroups=2, nGenes=5000, dropout=0.5, 
                     batchCells = c(400, 700), group.prob = c(0.3, 0.7),
                     de.prob = c(0.2, 0.1),
                     de.downProb = c(0.3, 0.4),
                     de.facLoc = 0.5,
                     de.facScale = 0.2){
  if (nGroups > 1) 
    method <- 'groups'
  else 
    method <- 'single'
  # group.prob <- rep(1, nGroups) / nGroups
  # new splatter requires dropout.type
  if ('dropout.type' %in% slotNames(newSplatParams())) 
  {
    if (dropout){
      dropout.type <- 'experiment'
    }  
    else{
      dropout.type <- 'none'
    }  
    sim <- splatSimulate(group.prob=group.prob, nGenes=nGenes, batchCells=batchCells,
                         dropout.type=dropout.type, method=method,seed=42, dropout.shape=-1, 
                         dropout.mid=dropout, 
                         de.prob=de.prob, de.downProb=de.downProb,
                         de.facLoc=de.facLoc, de.facScale=de.facScale
    )
  } else {
    sim <- splatSimulate(group.prob=group.prob, nGenes=nGenes, batchCells=batchCells,
                         dropout.present=!dropout, method=method,seed=42, dropout.shape=-1, 
                         dropout.mid=dropout,
                         de.prob=de.prob, de.downProb=de.downProb,
                         de.facLoc=de.facLoc, de.facScale=de.facScale
    ) 
  }
  print(class(sim))
  # splatter::splat
  #sim1 <- sim 
  #sim1@
  
  counts <- as.data.frame(t(counts(sim)))
  truecounts <- as.data.frame(t(assays(sim)$TrueCounts))
  dropout <- assays(sim)$Dropout
  mode(dropout) <- 'integer'
  dropout <- as.data.frame(t(dropout))
  cellinfo <- as.data.frame(colData(sim))
  geneinfo <- as.data.frame(rowData(sim))
  
  return(list(simu=sim,counts=counts,cellinfo=cellinfo,geneinfo=geneinfo,truecounts=truecounts,dropout=dropout))
}

sapply(1:dim(combi)[1],function(x){
  
  base_name <- paste0('simul',x,'_dropout_',gsub('\\.','',combi$d[[x]]),'_b1_',combi$b1[[x]],'_b2_',combi$b2[[x]],'/')
  dir.create(base_name, showWarnings = FALSE)
  
  sim <- simulate(dropout=combi$d[[x]], batchCells = c(combi$b1[[x]],combi$b2[[x]]))
  
  parameters_text<-file(paste0(base_name,"_parameters.txt"))
  writeLines(c(paste0("dropout = ",combi$d[[x]]),
               paste0("batchCells = c(",combi$b1[[x]],",",combi$b2[[x]],")"),
               "nGroups = 2",
               "nGenes = 5000",
               "group.prob = c(0.3, 0.7)",
               "de.prob = c(0.2, 0.1)",
               "de.downProb = c(0.3, 0.4)",
               "de.facLoc = 0.5",
               "de.facScale = 0.2"), parameters_text)
  close(parameters_text)

  counts <- sim$counts
  geneinfo <- sim$geneinfo
  cellinfo <- sim$cellinfo
  truecounts <- sim$truecounts
  dropout <- sim$dropout
  
  de_genes_ls <- rownames(geneinfo[(geneinfo$DEFacGroup1+geneinfo$DEFacGroup2)!=2,])
  de_genes_df <- geneinfo[de_genes_ls,]
  down_genes <- rownames(de_genes_df[de_genes_df$DEFacGroup1<de_genes_df$DEFacGroup2,])
  up_genes <- de_genes_ls[!de_genes_ls %in% down_genes]
  genes_deg <- c(up_genes, down_genes) 
  
  write.table(counts, file = paste0(base_name,"/counts.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(truecounts, file = paste0(base_name,"/truecounts.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(geneinfo, file = paste0(base_name,"/geneinfo.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(cellinfo, file = paste0(base_name,"/cellinfo.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(down_genes, file = paste0(base_name,"/true_down_genes.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(up_genes, file = paste0(base_name,"/true_up_genes.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  
})


# # Heatmap DEGs before batch effect removing 
# 
# counts_hm <- t(counts)
# colnames(counts_hm) <- paste(cellinfo[colnames(counts_hm),'Group'],colnames(counts_hm),sep="_")
# pbmc <- CreateSeuratObject(raw.data = counts_hm, project = "simulated_data", min.cells = 0, min.features = 200)
# pbmc
# pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e2)
# pbmc <- ScaleData(object = pbmc)
# 
# # 28/05/2019
# png("heatmap_before_batch_effect_removing/20190528_DEG_GT_all_hm.png", width = 2*480, height = 2*480, res = 2*72)
# gg <- DoHeatmap(object = pbmc, genes.use = genes_deg, slim.col.label = TRUE, group.label.rot =F,
#                remove.key = FALSE,col.low = "#0000FF",col.mid = "#FFFFFF", col.high = "#FF0000")
# # remove row labels 
# gg <- gg + theme(axis.text.y = element_blank())
# 
# # b <- ggplot_build(gg)
# # x_g1 <- median(seq(1:length(b$layout$panel_scales_x[[1]]$range$range)))
# # x_g2 <- median(seq(1:length(b$layout$panel_scales_x[[2]]$range$range)))
# x_g1 <- median(seq(1:length(pbmc@ident[pbmc@ident==levels(pbmc@ident)[1]])))
# x_g2 <- median(seq(1:length(pbmc@ident[pbmc@ident==levels(pbmc@ident)[2]])))
# y_g1 <- length(down_genes) + length(up_genes)/2
# y_g2 <- length(down_genes)/2
# 
# annotate_text <- data.frame(cell = c(x_g1,x_g2), 
#                    gene = c(y_g1,y_g2), 
#                    expression = c(0,0), 
#                    lab = c(length(up_genes),length(down_genes)),
#                    ident = c("Group1","Group2"))
# gg <- gg + geom_text(data=annotate_text, aes(x = cell, y = gene, label=lab),size=9)
# print(gg)
# dev.off()