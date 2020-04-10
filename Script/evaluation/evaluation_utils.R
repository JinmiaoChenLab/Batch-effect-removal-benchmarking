library(ggplot2)



###############################
# ARI Adjusted Rand Index
###############################
# Input: myData: 20 PCs and batch, celltype columns
# prefix_coln: the colnames of myData

ari_calcul_func <- function(myData, isOptimal=FALSE, 
                            method_use='resnet', prefix_coln='X_pca',
                            base_name='', 
                            maxcls=10, maxiter=30)
{
  library(cluster)
  library(mclust)
  set.seed(0)
  
  # Extract data
  pcs <- rep(1:20, 1)
  cpcs <- paste0(prefix_coln,pcs)
  # pca_mtx <- myData[,cpcs]
  
  
  cells_ls <- rownames(myData)
  ori_nbcells <- length(cells_ls)
  nbct <- length(unique(myData$celltype))
  percent_extract <- 0.8
  # Run func 10 times, each time extract 80% of data
  
  it <- c()
  total_ari_batch <- c()
  total_ari_batch_norm <- c()
  total_ari_celltype <- c()
  total_ari_celltype_norm <- c()
  total_fscoreARI <- c()
  nbiters <- 20
  for(i in 1:nbiters) {
    
    cells_extract <- sample(cells_ls, size=round(ori_nbcells*percent_extract), replace = F)
    
    myPCAExt <- myData[cells_extract,]
    
    ###############################
    # Clustering
    ###############################
    
    if(!isOptimal){  # isOptimal==FALSE
      
      # nbct : k equal number of unique cell types in the dataset
      clustering_result <- kmeans(x = myPCAExt[,cpcs], centers=nbct, iter.max = maxiter)
      
    } else if(isOptimal){
      #need to change kmeans default argument of max iteraction to >10, because
      #the dataset was shown to not converge with max iteration 10
      xtra_kmeans <- function(x, centers, nstart = 1L, 
                              algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", 
                                            "MacQueen"), trace = FALSE) {
        kmeans(x, centers, iter.max = maxiter, nstart = 1L, 
               algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"), 
               trace = FALSE)
      }
      # clusGap, run xtra_kmeans func first
      gapresult <- clusGap(x = myPCAExt[,cpcs], FUNcluster = xtra_kmeans,  
                           K.max = maxcls,
                           B = 100)
      optimal_k <- with(gapresult, maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
      clustering_result <- kmeans(x = myPCAExt[,cpcs], centers=optimal_k, 
                                  iter.max = maxiter)
    }
    
    myPCAExt$clusterlb <- clustering_result$cluster
    print('Nb clusters: ')
    print(length(unique(myPCAExt$clusterlb)))
    
    ###############################
    # ARI
    ###############################
    
    # run ARI
    ari_batch <- mclust::adjustedRandIndex(myPCAExt$batch, myPCAExt$clusterlb)
    ari_celltype<-mclust::adjustedRandIndex(myPCAExt$celltype, myPCAExt$clusterlb)
    
    # normalise values, ARI raw value from -1 to 1
    min_val <- -1
    max_val <- 1
    ari_batch_norm <- (ari_batch-min_val)/(max_val-min_val)
    ari_celltype_norm <- (ari_celltype-min_val)/(max_val-min_val)
    
    # produce final fscore ARI, similar to scMerge paper
    fscoreARI <- (2 * (1 - ari_batch_norm)*(ari_celltype_norm))/(1 - ari_batch_norm + ari_celltype_norm)
    
    it <- c(it,i)
    total_ari_batch <- c(total_ari_batch,ari_batch)
    total_ari_batch_norm <- c(total_ari_batch_norm,ari_batch_norm)
    total_ari_celltype <- c(total_ari_celltype,ari_celltype)
    total_ari_celltype_norm <- c(total_ari_celltype_norm, ari_celltype_norm)
    total_fscoreARI <- c(total_fscoreARI,fscoreARI)
    
  }  
  
  it <- c(it,nbiters+1)
  total_ari_batch <- c(total_ari_batch,round(median(total_ari_batch), digits = 3))
  total_ari_batch_norm <- c(total_ari_batch_norm,round(median(total_ari_batch_norm), digits = 3))
  total_ari_celltype <- c(total_ari_celltype,round(median(total_ari_celltype), digits = 3))
  total_ari_celltype_norm <- c(total_ari_celltype_norm, round(median(total_ari_celltype_norm), digits = 3))
  total_fscoreARI <- c(total_fscoreARI,round(median(total_fscoreARI), digits = 3))
  
  methods <- rep(method_use, nbiters)
  methods <- c(methods,paste0(method_use,'_median'))
  myARI <- data.frame("use_case"=methods, 
                      "iteration"=it,
                      "ari_batch"=total_ari_batch, 
                      "ari_batch_norm" = total_ari_batch_norm,
                      'ari_celltype' = total_ari_celltype,
                      'ari_celltype_norm' = total_ari_celltype_norm,
                      'fscoreARI' = total_fscoreARI)
  write.table(myARI, file = paste0(base_name,method_use,"_ARI.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  print('Save output in folder')
  print(base_name)
  
  return(myARI)
  
}  


# kBET evaluation
# input data frame: samples x features
kbet_calcul_func <- function(data, labelObs, base_name, label, sortlb, kn=20, flag_pca=FALSE)
{
  library(ggplot2)
  library(kBET)
  batch.estimate <- kBET(data, labelObs, k0=kn, plot=FALSE, do.pca = flag_pca)  
  
  plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                    each=length(batch.estimate$stats$kBET.observed)), 
                          data =  c(batch.estimate$stats$kBET.observed,
                                    batch.estimate$stats$kBET.expected))
  
  kbet_df <- data.frame("observed"=batch.estimate$stats$kBET.observed, "expected"=batch.estimate$stats$kBET.expected)
  # write.table(kbet_df,paste0(base_name,label,"_k",kn,"_evals.txt"), quote=F, sep='\t', row.names=FALSE, col.names=TRUE)
  
  ezTools::fast_save_table(kbet_df, base_name, paste0(label,"_k",kn,"_evals.txt"),
                           row.names = F, col.names = T, quote = F, sep = "\t")
  
  # write.table(plot.data,paste(base_name,label,"_k",kn,"_plot_data.txt",sep=""), quote=F, sep='\t', row.names=FALSE, col.names=TRUE)
  # save(batch.estimate, file = paste0(base_name,label,"_integrate_seurat_dataset1_uc3.rda"))
  # class(batch.estimate$stats)
  # load(paste0(base_name,"kbet_integrate_seurat_dataset1_uc3.rda"))
  
  g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
    labs(x='Test', y='Rejection rate',title=paste0('kBET test results ',sortlb)) +
    theme_bw() +  
    scale_y_continuous(limits=c(0,1))
  g
  library(ezTools)
  ezTools::save_image(g, paste0(base_name, label,"_k",kn,"bxplt.png"))
}  






# LISI_boxplot_fun(data = list(SeuratCCA_lisi_batch, Harmony_lisi_batch),vect_names_method = c('SeuratCCA','Harmony'), title='iLISI', base_name='', toptitle='Batch Mixing')
LISI_boxplot_fun <- function(data,vect_names_method, base_name, title='iLISI', toptitle='Batch Mixing', save_results=F){
  
  library(ggpubr)
  
  # combine the LISI index from all the method in one dataset
  if(length(data)>1){
    for (i in seq_along(data)){
      colnames(data[[i]]) <- c(vect_names_method[i],'cell')
    }
    data <- Reduce(function(x, y) merge(x, y, by="cell"), data)
    data <- data[,-1]
  } else {
    data <- as.data.frame(data[[1]][,1])
    # colnames(data) <- vect_names_method
  }
  
  cn <- c()
  for (i in seq_along(data)){
    cn <- c(cn, vect_names_method[i])
  }
  if(save_results){
    totaldf <- as.data.frame(data)
    colnames(totaldf) <- cn
    write.csv(totaldf, file=paste0(base_name, title,"_summary.csv"), row.names=FALSE)
    med_val <- c()
    for (c in cn){
      med_val <- c(med_val, round(median(totaldf[,c]), digits = 3))
    }
    median_df <- data.frame('methods'=cn, 'median_value'=med_val)
    write.csv(median_df, file=paste0(base_name, title,"_median_val.csv"), row.names=FALSE)
    print('Results saved')
  }
  
  # melt to plot 
  plotdata <- tidyr::gather(data)
  
  # put the y-axis in order 
  plotdata$key <- factor(plotdata$key,levels = rev(vect_names_method), order=T)
  # write.csv(plotdata, file = paste0(base_name,'summary_',title,'.csv'), row.names = T, col.names = T)
  # print('Save the output summary in the folder: ')
  # print(base_name)
  # plot
  
  p <- ggplot(plotdata, aes(key, value))
  p <- p + geom_boxplot(color='black',outlier.size = 1) + coord_flip() # outlier.size=-1 remove outliers dots
  p <- p + labs(title=toptitle,x=" ", y = title)
  p <- p + theme(axis.title.y = element_blank(), 
                 axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
                 axis.title.x = element_text(size=15),
                 axis.text.x = element_text(size=11,colour = 'black'),
                 axis.text.y = element_text(size=11,colour = 'black'),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(), 
                 plot.title = element_text(hjust = 0.5)) 
  # p <- p + theme_gray(base_size = 14)
  # axis.line.y=element_blank(), axis.ticks.y=element_blank()
  
  # save
  library(ezTools)
  ezTools::save_image(p, paste0(base_name, "boxplot_",title,".png"))
  print('Boxplot is saved')
  return(p)
}


scatter_plot_func <- function(plot_df, xstring, ystring, plottype, plottitle='Dataset', colorcode, xlabel='ARI(Cell type)', ylabel='1-ARI(batch)') {
  library(ggrepel)
  
  plot_df <- plot_df[order(as.character(plot_df[,plottype])),]
  
  plotobj <- ggplot(plot_df, aes_string(x = xstring, y = ystring, colour = plottype)) + 
    geom_point(shape=19, alpha = 1, size = 4) + theme_bw() #alpha = 0.6
  
  plotobj <- plotobj + geom_text_repel(aes_string(label = plottype),
                            nudge_y      = 0.02,
                            angle        = 0,
                            direction    = "y",
                            vjust        = 0,
                            segment.size = 0.2)
  #ggplot_build(plotobj)$data
  
  # xmax = max(plot_df[[xstring]])
  # xmin = min(plot_df[[xstring]])
  # ymax = max(plot_df[[ystring]])
  # ymin = min(plot_df[[ystring]])
  # xmax = ceiling(xmax)
  # xmin = floor(xmin)
  # ymax = ceiling(ymax)
  # ymin = floor(ymin)
  # 
  if(plottitle==''){
    plotobj <- plotobj + labs(x=xlabel,y=ylabel) 
  } else {
    plotobj <- plotobj + labs(x=xlabel,y=ylabel,title=plottitle) 
  }
  
  plotobj <- plotobj + theme(legend.title = element_blank(), 
                             plot.title = element_text(color="black", size=18, hjust = 0.5),
                             legend.position = "none", 
                             axis.line = element_line(color = "black", size = 0.5, linetype = "solid"), 
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.border = element_blank()
                             # axis.text = element_blank()
                             #axis.title = element_blank(),
                             # axis.ticks = element_blank()
  ) + scale_fill_manual(values = colorcode)
  
  # + scale_color_brewer(palette=marker_palette)
  
  return(plotobj)
}


runtime_plot <- function(totaldf, toptitle='Runtime', title='minutes'){
  
  prt <- ggplot(totaldf, aes(method_name, exetime_mins))
  prt <- prt + geom_bar(stat="identity", width=0.5,color="black", fill="white") + coord_flip()
  
  prt <- prt + labs(x=' ', y=title,title=toptitle)
  
  # prt <- prt + theme(axis.title.y = element_blank(),
  #                      axis.title.x = element_text(size=15),
  #                      axis.text.x = element_text(angle = 90,size=7,colour = 'black'),
  #                      axis.text.y = element_text(size=11,colour = 'black'),
  #                      plot.title = element_text(hjust = 0.5)) 
  prt <- prt + theme(axis.title.y = element_blank(), 
                     axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
                     axis.title.x = element_text(size=15),
                     axis.text.x = element_text(size=8,colour = 'black', angle = 90),
                     axis.text.y = element_text(size=11,colour = 'black'),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     panel.border = element_blank(), 
                     plot.title = element_text(hjust = 0.5)) 
  return(prt)
}




average_silhouette_width_plot <- function(asw_df, toptitle='Silhouette Width Batch Mixing', title='ASW'){
  
  pasw <- ggplot(asw_df, aes(method_name, asw))
  pasw <- pasw + geom_bar(stat="identity", width=0.5,color="black", fill="white") 
  pasw <- pasw + coord_flip()
  
  # p <- p + ylab(title) 
  
  pasw <- pasw + labs(x=' ', y=title,title=toptitle)
  
  pasw <- pasw + theme(axis.title.y = element_blank(), 
                       axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
                       axis.title.x = element_text(size=15),
                       axis.text.x = element_text(size=8,colour = 'black', angle = 90),
                       axis.text.y = element_text(size=11,colour = 'black'),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                       panel.border = element_blank(), 
                       plot.title = element_text(hjust = 0.5)) 
  # pasw <- pasw + theme(axis.title.y = element_blank(),
  #                      axis.title.x = element_text(size=15,family="TT Arial"),
  #                      axis.text.x = element_text(angle = 90,size=7,colour = 'black',family="TT Arial"),
  #                      axis.text.y = element_text(size=11,colour = 'black',family="TT Arial"),
  #                      plot.title = element_text(hjust = 0.5,family="TT Arial")) 
  return(pasw)
}  

# time util function
runtime_export_func <- function(t1, t2, labeluc, base_name, timeunits='secs_mins')
{
  label_usecase = labeluc  #label "usecase3_method_abc"
  time_secs = c()
  time_secs = c(time_secs, as.numeric(difftime(t2,t1,units='secs')))
  print(tail(time_secs,n=1))
  
  time_mins = c()
  time_mins = c(time_mins, as.numeric(difftime(t2,t1,units='mins')))
  print(tail(time_mins,n=1))
  
  myRunTime <- data.frame("use_case"=label_usecase, "exetime_secs"=time_secs, "exetime_mins"=time_mins)
  write.table(myRunTime, file = paste0(base_name,labeluc,"_runtime.txt"), row.names = FALSE, col.names = TRUE, sep="\t")
}  

# input: a data frame genes x cells
# base_name where to save the summary output 
# summary values list

summary_mtx <- function(normalized_df, base_name="", file_name="summary_normalized")
{
  normalized_df <- as.matrix(normalized_df)
  neg_percent <- length(which(normalized_df < 0)) * 100 / (nrow(normalized_df) * ncol(normalized_df))
  summary_vals <- c(min(normalized_df))
  summary_vals <- c(summary_vals, quantile(normalized_df, probs = c(.25, .5, .75)))
  summary_vals <- c(summary_vals, max(normalized_df), neg_percent)
  vals_name <- c('min','quantile25','quantile50','quantile75','max','percen_neg_vals')
  mySummary <- data.frame("stat"=vals_name, "vals"=summary_vals)
  write.table(mySummary, file = paste0(base_name,file_name,"_values.txt"), row.names = FALSE, col.names = TRUE, sep="\t")
  return(summary_vals)
}  


# min_val <- 0
# max_val <- nb batches or max_val = nb cell types

normalize_mtx <- function(myData,  cols_use='batch', min_val = 0, max_val=10, base_name='')
{
  med_val <- c()
  med_val_norm <- c()
  methods_use <- c()
  
  for(c in cols_use){
    methods_use <- c(methods_use, c)
    mi <- median(myData[,c])
    med_val <- c(med_val, mi)
    med_val_norm <- c(med_val_norm, (mi - min_val) / (max_val - min_val))
  }
  res <- data.frame('methods_use'=methods_use, 'iLISI_median'=med_val,
                            'iLISI_median_norm'=med_val_norm)
  
}  


lisi_perplexity_plot <- function(totaldf, xstring='method', 
                                 ystring='cLISI_median_celltype',
                                 plottype = 'perplexity_use',
                                 toptitle='Runtime', title='minutes'){
  
  prt <- ggplot(totaldf, aes_string(xstring, ystring, colour = plottype))
  prt <- prt + geom_line(stat="identity", width=0.5,color="black", fill="white") + coord_flip()
  prt <- prt + geom_point()
  prt <- prt + labs(x=' ', y=title,title=toptitle)
  
  # prt <- prt + theme(axis.title.y = element_blank(),
  #                      axis.title.x = element_text(size=15),
  #                      axis.text.x = element_text(angle = 90,size=7,colour = 'black'),
  #                      axis.text.y = element_text(size=11,colour = 'black'),
  #                      plot.title = element_text(hjust = 0.5)) 
  prt <- prt + theme(legend.title = element_text(size=11, hjust = 0.5),
                     axis.title.y = element_blank(), 
                     axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
                     axis.title.x = element_text(size=15),
                     axis.text.x = element_text(size=8,colour = 'black', angle = 90),
                     axis.text.y = element_text(size=11,colour = 'black'),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     panel.border = element_blank(), 
                     plot.title = element_text(hjust = 0.5)) 
  return(prt)
}


get_celltype_common_kbet_V2 <- function(data_dir, fn){
  
  
  if(file.exists(paste0(data_dir, fn))){
    
    myData <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
    colnames(myData)[grep('[cC]ell_?[tT]ype',colnames(myData))] <- 'cell_type'
    colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
    batches <- unique(myData$batch)
    celltypels <- unique(myData$cell_type)
    # print("List of cell type: ")
    # print(celltypels)
    if(length(batches)<=1){
      return(NULL)
    }
    
    ctls <- list()
    count <- 0
    for (b in batches){
      count <- count + 1
      ct <- unique(myData[which(myData$batch==b),'cell_type'])
      ctls[[count]] <- ct
      # print("Batch and cell type infos:")
      # print(b)
      # print(ct)
    }
    
    for(i in rep(1:length(ctls),1)){
      if(i==2){
        ct_common <- intersect(ctls[[i-1]], ctls[[i]]) 
      }
      if(i>2){   #more than 2 batches
        ct_common <- intersect(ct_common, ctls[[i]])    
      }
    }
    print('Common cells:')
    print(ct_common)
    cells_extract <- c()
    batchls <- c()
    for(ct in ct_common){
      print('ct :')
      print(ct)
      batch_cn <- unique(myData[which(myData$cell_type==ct),'batch'])
      batchls <- c(batchls, batch_cn)
      bls <- list()
      for(i in rep(1:length(batch_cn),1)){
        
        nbcells <- nrow(myData[which((myData$cell_type==ct) & (myData$batch==batch_cn[i])),])
        bls[[as.character(batch_cn[i])]] <- nbcells
        if(i==1){
          min_cells <- nbcells
        }
        if(min_cells > nbcells){
          min_cells <- nbcells
        }
      }
      # print("Min nb cells: ")
      # print(min_cells)
      
      for(b in batch_cn){
        print(b)
        print(bls[[as.character(b)]])
        cn <- rownames(myData[which((myData$cell_type==ct) & (myData$batch==b)),])
        if(bls[[as.character(b)]]>min_cells){
          set.seed(42)
          cn <- sample(cn, size=min_cells, replace = F)  
          
          # check the distribution of downsampling and whole population is similar
          
        }
        print(length(cn))
        cells_extract <- c(cells_extract, cn)
      } 
    }  
    print('nb extracted cells in total: ')
    print(length(cells_extract))
    print(summary(myData[cells_extract,'cell_type']))
    
    return(list('ct_common'=ct_common, 'cells_extract'=cells_extract, 
                'batches'=unique(batchls), 'ctls'=celltypels))
  }else{
    print("File does not exist, pls double check the directory")
    print(paste0(data_dir, fn))
  }  
  
}


get_celltype_common_kbet_ds <- function(data_dir, fn, percent_ext=0.1){
  
  
  if(file.exists(paste0(data_dir, fn))){
    
    myData <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
    colnames(myData)[grep('[cC]ell_?[tT]ype',colnames(myData))] <- 'cell_type'
    colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
    batches <- unique(myData$batch)
    celltypels <- unique(myData$cell_type)
    # print("List of cell type: ")
    # print(celltypels)
    if(length(batches)<=1){
      return(NULL)
    }
    
    ctls <- list()
    count <- 0
    for (b in batches){
      count <- count + 1
      ct <- unique(myData[which(myData$batch==b),'cell_type'])
      ctls[[count]] <- ct
      # print("Batch and cell type infos:")
      # print(b)
      # print(ct)
    }
    
    for(i in rep(1:length(ctls),1)){
      if(i==2){
        ct_common <- intersect(ctls[[i-1]], ctls[[i]]) 
      }
      if(i>2){   #more than 2 batches
        ct_common <- intersect(ct_common, ctls[[i]])    
      }
    }
    print('Common cells:')
    print(ct_common)
    cells_extract <- c()
    batchls <- c()
    for(ct in ct_common){
      print('ct :')
      print(ct)
      batch_cn <- unique(myData[which(myData$cell_type==ct),'batch'])
      batchls <- c(batchls, batch_cn)
      bls <- list()
      for(i in rep(1:length(batch_cn),1)){
        
        nbcells <- nrow(myData[which((myData$cell_type==ct) & (myData$batch==batch_cn[i])),])
        bls[[as.character(batch_cn[i])]] <- nbcells
        if(i==1){
          min_cells <- nbcells
        }
        if(min_cells > nbcells){
          min_cells <- nbcells
        }
      }
      # print("Min nb cells: ")
      # print(min_cells)
      
      for(b in batch_cn){
        print(b)
        print(bls[[as.character(b)]])
        cn <- rownames(myData[which((myData$cell_type==ct) & (myData$batch==b)),])
        if(bls[[as.character(b)]] > min_cells * percent_ext){
          set.seed(42)
          cn <- sample(cn, size=round(min_cells * percent_ext), replace = F)
          # cn <- sample(cn, size=min_cells, replace = F)  
          
          # check the distribution of downsampling and whole population is similar
          
        }
        print(length(cn))
        cells_extract <- c(cells_extract, cn)
      } 
    }  
    print('nb extracted cells in total: ')
    print(length(cells_extract))
    print(summary(myData[cells_extract,'cell_type']))
    
    return(list('ct_common'=ct_common, 'cells_extract'=cells_extract, 
                'batches'=unique(batchls), 'ctls'=celltypels))
  }else{
    print("File does not exist, pls double check the directory")
    print(paste0(data_dir, fn))
  }  
  
}

get_celltype_common_kbet <- function(data_dir, fn){
  myData <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
  colnames(myData)[grep('[cC]ell_?[tT]ype',colnames(myData))] <- 'cell_type'
  colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
  batches <- unique(myData$batch)
  celltypels <- unique(myData$cell_type)
  # print("List of cell type: ")
  # print(celltypels)
  if(length(batches)<=1){
    return(NULL)
  }
  
  ctls <- list()
  count <- 0
  for (b in batches){
    count <- count + 1
    ct <- unique(myData[which(myData$batch==b),'cell_type'])
    ctls[[count]] <- ct
    # print("Batch and cell type infos:")
    # print(b)
    # print(ct)
  }
  
  for(i in rep(1:length(ctls),1)){
    if(i==2){
      ct_common <- intersect(ctls[[i-1]], ctls[[i]]) 
    }
    if(i>2){   #more than 2 batches
      ct_common <- intersect(ct_common, ctls[[i]])    
    }
  }
  print('Common cells:')
  print(ct_common)
  cells_extract <- c()
  batchls <- c()
  for(ct in ct_common){
    print('ct :')
    print(ct)
    batch_cn <- unique(myData[which(myData$cell_type==ct),'batch'])
    batchls <- c(batchls, batch_cn)
    bls <- list()
    for(i in rep(1:length(batch_cn),1)){
      
      nbcells <- nrow(myData[which((myData$cell_type==ct) & (myData$batch==batch_cn[i])),])
      bls[[as.character(batch_cn[i])]] <- nbcells
      if(i==1){
        min_cells <- nbcells
      }
      if(min_cells > nbcells){
        min_cells <- nbcells
      }
    }
    # print("Min nb cells: ")
    # print(min_cells)
    
    for(b in batch_cn){
      print(b)
      print(bls[[as.character(b)]])
      cn <- rownames(myData[which((myData$cell_type==ct) & (myData$batch==b)),])
      if(bls[[as.character(b)]]>min_cells){
        set.seed(42)
        cn <- sample(cn, size=min_cells, replace = F)  
        
        # check the distribution of downsampling and whole population is similar
        
      }
      print(length(cn))
      cells_extract <- c(cells_extract, cn)
    } 
  }  
  print('nb extracted cells in total: ')
  print(length(cells_extract))
  print(summary(myData[cells_extract,'cell_type']))
  
  return(list('ct_common'=ct_common, 'cells_extract'=cells_extract, 
              'batches'=unique(batchls), 'ctls'=celltypels))
  
}


run_kBET_final <- function(cells_extract, fn, data_dir, save_dir, eval_metric, methods_use, kn=30){
  if(file.exists(paste0(data_dir, fn))){
    myData <- ezTools::fast_read_table(paste0(data_dir, fn), row.names = T, 
                                       col.names = T, sep = ",", check.names = FALSE)
    
    s = rownames(myData) %in% cells_extract
    if(sum(s==T)>0){
      myData <- myData[s,]
      print(paste0('Nb extracted cells: ', nrow(myData), ' ', methods_use))
      cpcs <- grep('([Pp][Cc]_?)|(V)|(Harmony)|(W)',colnames(myData))
      cpcs <- cpcs[1:20]
      
      colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
      
      dir.create(paste0(save_dir, eval_metric, methods_use,'/'),showWarnings = FALSE)
      print(cat('Nb batches: ', unique(myData$batch),"\n",sep="\t"))
      kbet_calcul_func(data=myData[,cpcs], labelObs=myData$batch, 
                       base_name=paste0(save_dir, eval_metric, methods_use,'/'), 
                       label="kbet_batch",sortlb="batch", kn)  
    }else{
      print("****************************")
      print('Cell name error, check again')
    }
    
  }else{
    print("************* BUG ***************")
    print("File does not exist, double check the directory again")
    print(methods_use)
  }
}  




summary_KBET <- function(data_dir, methods_use, method_dir, dataset_use, eval_metric, save_dir, kn=30){
  fn <- paste0(dataset_use, '_raw_pca.csv')
  if(file.exists(paste0(data_dir, fn))){
    meta_ls <- get_celltype_common_kbet(data_dir, fn)
    if(length(meta_ls$ct_common)>0){
      
      fn_ls <- c()
      for(i in rep(1:length(method_dir),1)){
        fn_ls <- c(fn_ls, paste0(dataset_use,'_', method_dir[i],'_pca','.csv'))
      }
      print(fn_ls)
      
      for(i in rep(1:length(methods_use), 1)){
        print(methods_use[i])
        run_kBET_final(meta_ls$cells_extract, fn_ls[i], data_dir, this_dir, eval_metric, methods_use[i], kn)
      }
    }
  }
  return(meta_ls)
} 


summary_KBET_v2 <- function(meta_ls, data_dir, methods_use, method_dir, dataset_use, eval_metric, save_dir, kn=30){
  if(length(meta_ls$ct_common)>0){
    
    fn_ls <- c()
    for(i in rep(1:length(method_dir),1)){
      fn_ls <- c(fn_ls, paste0(dataset_use,'_', method_dir[i],'_pca','.csv'))
    }
    print(fn_ls)
    
    for(i in rep(1:length(methods_use), 1)){
      print(methods_use[i])
      run_kBET_final(meta_ls$cells_extract, fn_ls[i], data_dir, this_dir, eval_metric, methods_use[i], kn)
    }
  }
} 


summary_KBET_d9 <- function(meta_ls, data_dir, methods_use, method_dir, dataset_use, eval_metric, save_dir, kn=30){
  
  fn_ls <- c()
  for(i in rep(1:length(method_dir),1)){
    fn_ls <- c(fn_ls, paste0(dataset_use,'_', method_dir[i],'_pca','.csv'))
  }
  print(fn_ls)
  
  for(i in rep(1:length(methods_use), 1)){
    print(methods_use[i])
    run_kBET_final(meta_ls$cells_extract, fn_ls[i], data_dir, this_dir, eval_metric, methods_use[i], kn)
  }
  
} 


