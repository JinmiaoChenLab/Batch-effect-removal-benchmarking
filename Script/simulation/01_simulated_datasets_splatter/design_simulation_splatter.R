
# Tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splat_params.R

rm(list=ls())

library(splatter)

this.dir <- "../../../Data/dataset3/"
setwd(this.dir)

# Ex: ./simul1_dropout_005_b1_500_b2_900/    unbalanced number of cells in 2 batches, small dropout
# Ex: ./simul2_dropout_025_b1_500_b2_900/    unbalanced number of cells in 2 batches, large dropout
# Ex: ./simul3_dropout_005_b1_500_b2_450/    balanced number of cells in 2 batches, small dropout
# Ex: ./simul3_dropout_025_b1_500_b2_450/    unbalanced number of cells in 2 batches, large dropout
# Ex: ./simul4_dropout_005_b1_80_b2_400/     small number of cells in 2 batches, small dropout
# Ex: ./simul4_dropout_025_b1_80_b2_400/     small number of cells in 2 batches, large dropout

dropout <- rep(c(0.05,0.25),times=3)
b1 <- rep(c(500,500,80),each = 2)
b2 <- rep(c(900,450,400),each = 2)
combi <- data.frame(d=dropout,b1=b1,b2=b2)

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
