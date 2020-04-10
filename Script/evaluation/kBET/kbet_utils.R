normalize_values <- function(myData, colns, min_val=1, max_val=2){
  med_val <- c()
  med_val_norm <- c()
  methods_use <- c()
  for(c in colnames(myData)){
    methods_use <- c(methods_use, c)
    mi <- median(myData[,c])
    med_val <- c(med_val, mi)
    med_val_norm <- c(med_val_norm, (mi - min_val) / (max_val - min_val))
  }
  myDataNorm <- data.frame('X1'=methods_use, 'X2'=med_val, 'X3'=med_val_norm)
  colnames(myDataNorm) <- colns
  return(myDataNorm)
}


get_celltype_common <- function(data_dir, fn){
  myData <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
  colnames(myData)[grep('[cC]ell_?[tT]ype',colnames(myData))] <- 'cell_type'
  colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
  print(dim(myData))
  batches <- unique(myData$batch)
  celltypels <- unique(myData$cell_type)
  print(celltypels)
  if(length(batches)>1){
    ctls <- list()
    count <- 0
    for (b in batches){
      count <- count + 1
      ct <- unique(myData[which(myData$batch==b),'cell_type'])
      ctls[[count]] <- ct
      print("Batch")
      print(b)
      print(ct)
    }
    for(i in rep(1:length(ctls))){
      if(i==2){
        ct_common <- intersect(ctls[[i-1]], ctls[[i]])    
      }
      if(i>2){   #more than 2 batches
        ct_common <- intersect(ct_common, ctls[[i]])    
      }
    }
    ct_common <- unique(ct_common)
    print(ct_common)
    cells_common <- rownames(myData)[which(myData$cell_type %in% ct_common) ]
    return(list('ct_common'=ct_common, 'cells_common'=cells_common, 'batches'=batches, 'celltypels'=celltypels))
  } else{
    return(NULL)
  } 
  
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




#####################################
# kBET collect all results
#####################################
get_kbet_output <- function(dataset_use, save_dir='', kn = 30, eval_metric = 'kBET/',save_results=T){
  
  seurat2_df <- read.table(paste0(eval_metric,"Seurat_2/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  seurat3_df <- read.table(paste0(eval_metric,"Seurat_3/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  harmony_df <- read.table(paste0(eval_metric,"Harmony/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  fastMNN_df <- read.table(paste0(eval_metric,"fastMNN/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  resnet_df <- read.table(paste0(eval_metric,"MMD-ResNet/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  scanorama_df <- read.table(paste0(eval_metric,"Scanorama/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  scGen_df <- read.table(paste0(eval_metric,"scGen/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)  
  raw_df <- read.table(paste0(eval_metric,"Raw/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)  
  correctMNN_df <- read.table(paste0(eval_metric,"MNN_Correct/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)  
  combat_df <- read.table(paste0(eval_metric,"ComBat/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE) 
  liger_df <- read.table(paste0(eval_metric,"LIGER/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE) 
  limma_df <- read.table(paste0(eval_metric,"limma/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE) 
  scMerge_df <- read.table(paste0(eval_metric,"scMerge/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE) 
  zinbwave_df <- read.table(paste0(eval_metric,"ZINB-WaVE/kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE) 
  
  
  methods_use <- c('Raw','Seurat_2','Seurat_3','Harmony','fastMNN','MNN_Correct','ComBat',
                   'limma','scGen','Scanorama','MMD-ResNet','ZINB-WaVE','scMerge','LIGER')
  dir.create(paste0(save_dir, eval_metric,"result/"), showWarnings = F)
  dir.create(paste0(save_dir, eval_metric,"result/",kn,"/"), showWarnings = F)
  pkbet <- kbet_plot_func(data = list(raw_df, seurat2_df, seurat3_df,
                                         harmony_df, fastMNN_df, correctMNN_df,
                                         combat_df, limma_df, scGen_df,                    
                                         scanorama_df, resnet_df, 
                                         zinbwave_df, scMerge_df, liger_df),
                             vect_names_method = methods_use, 
                             title='Acceptance rate', toptitle=paste0('kBET k ',kn, ' ',dataset_use),
                             base_name=paste0(save_dir, eval_metric, "result/",kn,"/"),save_results=T)
  
  
  
  save(pkbet, file = paste0(save_dir, eval_metric, "result/",kn, "/pkbet_",dataset_use,"_plots_",kn,".rda"))
  return(pkbet)
}



kbet_plot_func <- function(data, vect_names_method, base_name, 
                           title='Rejection rate',toptitle='kBET batch effect test', 
                           save_results=F)
{
  for (i in seq_along(data))
  {
    if(i==1)
    {
      totaldf <- data[[1]]$observed
      total_acceptance <- 1 - data[[1]]$observed
    }  
    else{
      totaldf <- cbind(totaldf, data[[i]]$observed)
      total_acceptance <- cbind(total_acceptance, 1 - data[[i]]$observed)
    }
  }
  
  cn <- c()
  for (i in seq_along(data)){
    cn <- c(cn, vect_names_method[i])
  }
  totaldf <- as.data.frame(totaldf)
  total_acceptance <- as.data.frame(total_acceptance)
  colnames(totaldf) <- cn
  colnames(total_acceptance) <- cn
  if(save_results){
    write.csv(totaldf, file=paste0(base_name,"total_kBET_observed.csv"), row.names=FALSE)
    write.csv(total_acceptance, file=paste0(base_name,"total_kBET_observed_acceptance.csv"), row.names=FALSE)
    med_val <- c()
    for (c in cn){
      med_val <- c(med_val, round(median(totaldf[,c]), digits = 3))
    }
    median_df <- data.frame('methods'=cn, 'median_observed'=med_val, 'median_acceptance_rate'=(1-med_val))
    write.csv(median_df, file=paste0(base_name, title,"_median_val_kBET.csv"), row.names=FALSE)
  }  
  
  # View(totaldf)
  plotdata <- tidyr::gather(total_acceptance)
  
  # put the y-axis in order 
  plotdata$key <- factor(plotdata$key,levels = rev(vect_names_method), order=T)
  # write.csv(plotdata, file = paste0(base_name,'summary_kBET.csv'), row.names = T, col.names = T)
  # print('Save the output summary in the folder: ')
  # print(base_name)
  
  p <- ggplot(plotdata, aes(key, value))
  p <- p + geom_boxplot(color='black') + coord_flip()
  # p <- p + ylab(title) 
  p <- p + labs(x=' ', y=title,title=toptitle)
  p <- p + theme(axis.title.y = element_blank(), 
                 axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
                 axis.title.x = element_text(size=15),
                 axis.text.x = element_text(size=14,colour = 'black'),
                 axis.text.y = element_text(size=11,colour = 'black'),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(), 
                 plot.title = element_text(hjust = 0.5)) 
  # p <- p + theme(axis.title.y = element_blank(),
  #                axis.title.x = element_text(size=15,family="TT Arial"),
  #                axis.text.x = element_text(size=12,colour = 'black',angle = 90,family="TT Arial"),
  #                axis.text.y = element_text(size=11,colour = 'black',family="TT Arial"),
  #                plot.title = element_text(hjust = 0.5,family="TT Arial")) 
  
  # p <- p + theme_gray(base_size = 14)
  # supprimes les outliers en utilisant le paramètre outlier.size=-1 dans geom_boxplot() et modifier les limites de l’axe x avec xlim(). 
  # mets des points plus petits en utilisant le parameter outlier.size = 5 dans geom_boxplot() 
  # axis.line.y=element_blank(), axis.ticks.y=element_blank(),
  # p
  # save
  # p
  library(ezTools)
  ezTools::save_image(p, paste0(base_name, "boxplot_",title,".png"))
  print('Boxplot is saved')
  return(p)
}  

generate_median_kbet <- function(pct, kns, res_dir, dataset_use, data_dir='',kbet_output=T){
  uc <- paste0(pct,"%")
  for(i in 1:length(kns)){
    kn = kns[i]
    print(paste0("k0 values:  ", kn))
    # Get all median values of kbet for each k and all methods
    if(kbet_output){
      pkbet <- get_kbet_output(dataset_use, data_dir, kn, eval_metric, save_results=T)
    }
    
    dir_result <- paste0(res_dir, kn, "/Acceptance rate_median_val_kBET.csv")
    
    if(file.exists(dir_result)){
      res_i = read.csv(dir_result, header=T)
      print(colnames(res_i))
      rej_i <- subset(res_i, select = c("methods","median_observed"))
      acc_i <- subset(res_i, select = c("methods","median_acceptance_rate"))
      colnames(rej_i)[which(colnames(rej_i) == 'median_observed')] <- paste0(pct[i],"%")
      colnames(acc_i)[which(colnames(acc_i) == 'median_acceptance_rate')] <- paste0(pct[i],"%")
      if(i==1){
        rej = rej_i
        acc = acc_i
      }else{
        rej = merge(rej, rej_i, by="methods")  
        acc = merge(acc, acc_i, by="methods") 
      }
    }  
  }
  
  # colnames(rej) <- c("methods", uc)
  # colnames(acc) <- c("methods", uc)
  med_rej <- c()
  med_acc <- c()
  for(r in 1:nrow(rej)){
    med_rej <- c(med_rej, median(as.numeric(rej[r,uc])))
    med_acc <- c(med_acc, median(as.numeric(acc[r,uc])))
  }
  rej$median_val <- med_rej
  acc$median_val <- med_acc
  
  rej = rej[order(rej$median_val, decreasing = F),]
  acc = acc[order(acc$median_val, decreasing = T),]
  
  write.csv(rej, paste0(res_dir,dataset_use,"_kbet_rejection.csv"))
  write.csv(acc,paste0(res_dir,dataset_use,"_kbet_acception.csv"))
  
}


get_equal_nbcells_d9 <- function(data_dir, fn){
  if(file.exists(paste0(data_dir, fn))){
    myData <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
    colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
    batches <- unique(myData$batch)
    if(length(batches)<=1){
      return(NULL)
    }
    bls <- list()
    cells_extract <- c()
    for(i in rep(1:length(batches),1)){
      nbcells <- nrow(myData[which(myData$batch==batches[i]),])
      bls[[as.character(batches[i])]] <- nbcells
      if(i==1){
        min_cells <- nbcells
      }
      if(min_cells > nbcells){
        min_cells <- nbcells
      }
    }
    for(b in batches){
      print('_______________')
      print(b)
      print(bls[[as.character(b)]])
      cn <- rownames(myData[which(myData$batch==b),])
      if(bls[[as.character(b)]]>min_cells){
        cn <- sample(cn, size=min_cells, replace = F)  
        # check the distribution of downsampling and whole population is similar
      }
      print(length(cn))
      cells_extract <- c(cells_extract, cn)
    } 
    print('Nb extracted cells in total: ')
    print(length(cells_extract))
    print(summary(myData[cells_extract,'batch']))
    return(list('cells_extract'=cells_extract, 'batches'=unique(batches)))
  }else{
    print("File does not exist, pls double check the directory")
    print(paste0(data_dir, fn))
  }  
}  


summary_KBET_d6 <- function(meta_ls_13, meta_ls_23, data_dir, methods_use, method_dir, dataset_use, eval_metric, save_dir, kn){
  
  
  print(paste0('K0 value use: ', kn))
  fn_ls <- c()
  for(i in rep(1:length(method_dir),1)){
    fn_ls <- c(fn_ls, paste0(dataset_use,'_', method_dir[i],'_pca','.csv'))
  }
  print(fn_ls)
  
 
  if(length(meta_ls_13$ct_common)>0){
    for(i in rep(1:length(methods_use), 1)){
      print(methods_use[i])
      run_kBET_final_d6('13', meta_ls_13$cells_extract, fn_ls[i], 
                        data_dir, this_dir, eval_metric, methods_use[i], kn)
    }
  }
  
  if(length(meta_ls_23$ct_common)>0){
    for(i in rep(1:length(methods_use), 1)){
      print(methods_use[i])
      run_kBET_final_d6('23', meta_ls_13$cells_extract, fn_ls[i], 
                        data_dir, this_dir, eval_metric, methods_use[i], kn)
    }
  }
  
  get_kbet_output_d6('13', dataset_use, save_dir, kn, eval_metric, save_results=T)
  get_kbet_output_d6('23', dataset_use, save_dir, kn, eval_metric, save_results=T)
  
  
  my13kbet <- read.csv(paste0(save_dir, eval_metric,"result/",kn,"/",'13','/','total_kBET_observed_acceptance.csv'), 
                       head=T, check.names = FALSE)
  
  my23kbet <- read.csv(paste0(save_dir, eval_metric,"result/",kn,"/",'23','/','total_kBET_observed_acceptance.csv'), 
                       head=T, check.names = FALSE)
  ls_kbet <- list()
  for(c in colnames(my13kbet)){
    ls_kbet[[c]] <- (my13kbet[,c] + my23kbet[,c]) / 2
  }
  mytotal_kbet <- as.data.frame(ls_kbet)
  
  # mytotal_kbet <- rbind(my13kbet, my23kbet)
  write.csv(mytotal_kbet, file=paste0(save_dir, eval_metric,"result/",kn,"/","total_kBET_observed_acceptance.csv"), row.names=FALSE)
  med_val <- c()
  cn <- colnames(mytotal_kbet)
  for (c in cn){
    med_val <- c(med_val, round(median(mytotal_kbet[,c]), digits = 3))
  }
  median_df <- data.frame('methods'=cn, 'median_acceptance_rate'=med_val,'median_observed'=(1-med_val))
  median_df <- median_df[order(median_df$median_acceptance_rate, decreasing = T),]
  write.csv(median_df, file=paste0(save_dir, eval_metric,"result/",kn,"/", "kbet_median_acceptance_rate_total.csv"), row.names=FALSE)
} 


get_common_celltype_d6 <- function(observed_bls, data_dir='', fn='raw_pca.csv'){
  if(file.exists(paste0(data_dir, fn))){
    myData <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
    colnames(myData)[grep('[cC]ell_?[tT]ype',colnames(myData))] <- 'cell_type'
    colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
    batches <- unique(myData$batch)
    celltypels <- unique(myData$cell_type)
    
    if(sum(observed_bls %in% batches == TRUE)==length(observed_bls)){
      batchls <- observed_bls
    }else if(sum(as.character(observed_bls) %in% batches == TRUE)==length(observed_bls)){
      batchls <- as.character(observed_bls)
    }else{
      print('Check again input of batches label')
    }
    if(length(batchls)>1){
      ctls <- list()
      count <- 0
      for (b in batchls){
        count <- count + 1
        ct <- unique(myData[which(myData$batch==b),'cell_type'])
        ctls[[count]] <- ct
      }
      for(i in rep(1:length(ctls))){
        if(i==2){
          ct_common <- intersect(ctls[[i-1]], ctls[[i]])    
        }
        if(i>2){   #more than 2 batches
          ct_common <- intersect(ct_common, ctls[[i]])    
        }
      }
    }  
    print(paste0('Observed batches: ', observed_bls))
    print(paste0('Cell type in common: ', ct_common))
    cells_extract <- c()
    batchesls <- c()
    for(ct in ct_common){
      print('ct :')
      print(ct)
      batch_cn <- unique(myData[which(myData$cell_type==ct),'batch'])
      batchesls <- c(batchesls, batch_cn)
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
                'batches'=unique(batchesls), 'ctls'=celltypels))
    
  }
  
}

run_kBET_final_d6 <- function(observed_blb, cells_extract, fn, data_dir, save_dir, eval_metric, methods_use, kn=30){
  if(file.exists(paste0(data_dir, fn))){
    myData <- ezTools::fast_read_table(paste0(data_dir, fn), row.names = T, 
                                       col.names = T, sep = ",", check.names = FALSE)
    
    cells_use <- rownames(myData) %in% cells_extract
    if(sum(cells_use==T) > 0){
      myData <- myData[cells_use,]
      print(dim(myData))
      cpcs <- grep('([Pp][Cc]_?)|(V)|(Harmony)',colnames(myData))
      cpcs <- cpcs[1:20]
      
      colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
      
      dir.create(paste0(save_dir, eval_metric, methods_use,'/'),showWarnings = FALSE)
      
      kbet_calcul_func(data=myData[,cpcs], labelObs=myData$batch, 
                       base_name=paste0(save_dir, eval_metric, methods_use,'/'), 
                       label=paste0(observed_blb,"_kbet_batch"),sortlb="batch", kn)  
    }else{
      print("****************************")
      print('Cell name error, check again')
    }
    
  }else{
    print("File does not exist, double check the directory again")
    print(methods_use)
  }
}  

get_downsample_d9 <- function(data_dir, fn, percent_ext=0.1){
  if(file.exists(paste0(data_dir, fn))){
    myData <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
    colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
    batches <- unique(myData$batch)
    if(length(batches)<=1){
      return(NULL)
    }
    bls <- list()
    cells_extract <- c()
    for(i in rep(1:length(batches),1)){
      nbcells <- nrow(myData[which(myData$batch==batches[i]),])
      bls[[as.character(batches[i])]] <- nbcells
      if(i==1){
        min_cells <- nbcells
      }
      if(min_cells > nbcells){
        min_cells <- nbcells
      }
    }
    for(b in batches){
      print('_______________')
      print(b)
      print(bls[[as.character(b)]])
      cn <- rownames(myData[which(myData$batch==b),])
      if(bls[[as.character(b)]] > min_cells * percent_ext){
        cn <- sample(cn, size=round(min_cells * percent_ext), replace = F)  
        # check the distribution of downsampling and whole population is similar
      }
      print(length(cn))
      cells_extract <- c(cells_extract, cn)
    } 
    print('Nb extracted cells in total: ')
    print(length(cells_extract))
    print(summary(myData[cells_extract,'batch']))
    return(list('cells_extract'=cells_extract, 'batches'=unique(batches)))
  }else{
    print("File does not exist, pls double check the directory")
    print(paste0(data_dir, fn))
  }  
}  
#####################################
# kBET collect all results
#####################################
get_kbet_output_d6 <- function(blabel, dataset_use, save_dir='', kn = 30, eval_metric = 'kBET/',save_results=T){
  
  seurat2_df <- read.table(paste0(eval_metric,"Seurat_2/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  seurat3_df <- read.table(paste0(eval_metric,"Seurat_3/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  harmony_df <- read.table(paste0(eval_metric,"Harmony/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  fastMNN_df <- read.table(paste0(eval_metric,"fastMNN/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  resnet_df <- read.table(paste0(eval_metric,"MMD-ResNet/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  scanorama_df <- read.table(paste0(eval_metric,"Scanorama/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)
  scGen_df <- read.table(paste0(eval_metric,"scGen/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)  
  raw_df <- read.table(paste0(eval_metric,"Raw/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)  
  correctMNN_df <- read.table(paste0(eval_metric,"MNN_Correct/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE)  
  combat_df <- read.table(paste0(eval_metric,"ComBat/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE) 
  liger_df <- read.table(paste0(eval_metric,"LIGER/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE) 
  limma_df <- read.table(paste0(eval_metric,"limma/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE) 
  scMerge_df <- read.table(paste0(eval_metric,"scMerge/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE) 
  zinbwave_df <- read.table(paste0(eval_metric,"ZINB-WaVE/",blabel,"_kbet_batch_k",kn,"_evals.txt"), head=T, check.names = FALSE) 
  
  
  methods_use <- c('Raw','Seurat_2','Seurat_3','Harmony','fastMNN','MNN_Correct','ComBat',
                   'limma','scGen','Scanorama','MMD-ResNet','ZINB-WaVE','scMerge','LIGER')
  dir.create(paste0(save_dir, eval_metric,"result/"), showWarnings = F)
  dir.create(paste0(save_dir, eval_metric,"result/",kn,"/"), showWarnings = F)
  dir.create(paste0(save_dir, eval_metric,"result/",kn,"/",blabel,'/'), showWarnings = F)
  pkbet <- kbet_plot_func(data = list(raw_df, seurat2_df, seurat3_df,
                                      harmony_df, fastMNN_df, correctMNN_df,
                                      combat_df, limma_df, scGen_df,                    
                                      scanorama_df, resnet_df, 
                                      zinbwave_df, scMerge_df, liger_df),
                          vect_names_method = methods_use, 
                          title='Acceptance rate', toptitle=paste0('kBET k ',kn, ' ',dataset_use,' b', blabel),
                          base_name=paste0(save_dir, eval_metric, "result/",kn,"/",blabel,'/'),save_results=T)
  
  
  
  save(pkbet, file = paste0(save_dir, eval_metric, "result/",kn, "/pkbet_",dataset_use,'_',blabel,"_plots_",kn,".rda"))
  return(pkbet)
}



######################## 
# other dataset: intersection
# d6: special, use union
######################## 

get_all_celltype_d6 <- function(data_dir, fn){
  myData <- read.csv(paste0(data_dir,'/',fn), head=T, row.names = 1, check.names = FALSE)
  colnames(myData)[grep('[cC]ell_?[tT]ype',colnames(myData))] <- 'cell_type'
  colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
  batches <- unique(myData$batch)
  if(length(batches)>1){
    ctls <- list()
    count <- 0
    for (b in batches){
      count <- count + 1
      ct <- unique(myData[which(myData$batch==b),'cell_type'])
      ctls[[count]] <- ct
    }
    for(i in rep(1:length(ctls))){
      if(i==2){
        ct_common <- union(ctls[[i-1]], ctls[[i]])    
      }
      if(i>2){   #more than 2 batches
        ct_common <- union(ct_common, ctls[[i]])    
      }
    }
    ct_common <- unique(ct_common)
    cells_common <- rownames(myData)[which((myData$batch %in% batches)
                                           & (myData$cell_type %in% ct_common)) ]
    return(list('ct_common'=ct_common, 'cells_common'=cells_common, 'batches'=batches))
  } else{
    return(NULL)
  } 
  
}
