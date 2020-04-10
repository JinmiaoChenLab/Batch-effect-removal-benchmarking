library(reshape2)
library(xlsx)

rm(list=ls())
setwd('/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/')

utils_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/source_code_github/evaluation_utils/'
source(paste0(utils_dir,'evaluation_utils.R'))
source(paste0(utils_dir,'lisi/lisi_utils.R'))


# Author : Marion Chevrier
# Date : 27/08/2019
# Aim : Statistical test to compare the evaluation methods (kBET, LISI, ARI, ASW) of the 13 methods (without BBKNN) + raw on all datasets 
# Updated: Hoa Tran


dataset_vect <- paste0('dataset',c(1:10)[-3])
kbet_stat_test <- TRUE
lisi_stat_test <- TRUE
ari_stat_test <- TRUE
asw_stat_test <- TRUE


sapply(dataset_vect, function(dataset){
  
  print(dataset)
  save_dir <- '/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/manuscript_results/evaluation/lisi_Hoa_v2/'
  dataset2 <- gsub('dataset','d',dataset)
  
  dir.create(paste0(save_dir, dataset2,'/statistical_test/'),showWarnings = FALSE)
  
  
  # obtain kBET for 5, 10, 15, 20 and 25%
  dataset_g <- ifelse(dataset=='dataset1','dataset1_uc3',
                      ifelse(dataset=='dataset2','dataset2_cellatlas',
                             ifelse(dataset=='dataset4','dataset4_human_pancreatic',
                                    ifelse(dataset=='dataset5','dataset5_human_pbmc',
                                           ifelse(dataset=='dataset6','dataset6_cell_line',
                                                  ifelse(dataset=='dataset7','dataset7_Mouse_retina',
                                                         ifelse(dataset=='dataset8','dataset8_Mouse_brain',
                                                                ifelse(dataset=='dataset9','dataset9_Human_cell_atlas',
                                                                       ifelse(dataset=='dataset10','dataset10_hematoMNN_Hoa',NA)))))))))

  dataset_h <- paste0('generate_PCA_tSNE_UMAP/',dataset,'/')

  this_dir <- paste0('/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/manuscript_results/',dataset_g,'/')
  data_dir = paste0('/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/',dataset_h)
  
  ################
  ## KBET TEST
  ################
  if(kbet_stat_test){
    print('kBET')
    fn <- paste0(dataset, '_raw_pca.csv')
    folder_kBEt <- 'kBET_Hoa_v2/'
    result_file <- '/total_kBET_observed_acceptance.csv'
    if(dataset=='dataset6'){
      meta_ls_13 <- get_common_celltype_d6(c(1,3), data_dir, fn)
      meta_ls_23 <- get_common_celltype_d6(c(2,3), data_dir, fn)
      size = (length(meta_ls_13$cells_extract) + length(meta_ls_23$cells_extract))/2

    } else if(dataset=='dataset9'){
      meta_ls <- get_downsample_d9(data_dir, fn, percent_ext=0.2)
      size = length(meta_ls$cells_extract)

    } else if((dataset=='dataset7') | (dataset=='dataset8')){
      meta_ls <- get_celltype_common_kbet_ds(data_dir, fn, percent_ext=0.2)
      size = length(meta_ls$cells_extract)
    } else {
      meta_ls <- get_celltype_common_kbet_V2(data_dir, fn)
      size = length(meta_ls$cells_extract)
    }

    pct = c(seq(from = 5, to = 25, by = 5))
    kns = floor(pct*size/100)

    res_dir = paste0(this_dir, folder_kBEt,'result/')

    kbet_result_list <- lapply(kns,function(l){
      read.csv(paste0(res_dir,l,result_file), check.names = F)
    })
    kbet_result <- do.call(rbind, kbet_result_list)
    kbet_result <- melt(kbet_result)

    # the data are not normally distributed
    #shapiro.test(kbet_result$value[kbet_result$variable == "Harmony"])# p-value < 2.2e-16
    #hist(kbet_result$value[kbet_result$variable == "Harmony"])
    # so prefer to use wilcoxon test instead of t-test
    test_wilcoxon <- pairwise.wilcox.test(kbet_result$value, kbet_result$variable, p.adj = "BH")

    table_pvalue <- test_wilcoxon$p.value
    table_pvalue_round <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'<0.001',round(table_pvalue,3)))
    table_pvalue_round <- data.frame(table_pvalue_round, check.names = FALSE)
    write.xlsx(x = table_pvalue_round, file = paste0(save_dir, dataset2,'/statistical_test/kBET_wilcoxon_test.xlsx'),showNA = F,sheetName = 'pvalue')
    table_pvalue2 <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'***',ifelse(table_pvalue<0.01,'**',ifelse(table_pvalue<0.05,'*','NS'))))
    table_pvalue2 <- data.frame(table_pvalue2, check.names = FALSE) 
    write.xlsx(x = table_pvalue2,file = paste0(save_dir, dataset2,'/statistical_test/kBET_wilcoxon_test.xlsx'),showNA = F,sheetName = 'pvalue2', append=TRUE)
  }
  # ###################
  # # LISI TEST
  # ###################
  if(lisi_stat_test){
    print('iLISI')
    ilisi_result <- read.csv(paste0(save_dir, dataset2,'/iLISI_summary.csv'), check.names = F)
    # if(dataset=='dataset6'){  # have updated, sep=','
    #   ilisi_result <- read.csv(paste0(save_dir, dataset2,'/iLISI_summary.csv'), sep='\t')
    # } else {
    #   ilisi_result <- read.csv(paste0(save_dir, dataset2,'/iLISI_summary.csv'))
    # }
    ilisi_result <- melt(ilisi_result)
    test_wilcoxon <- pairwise.wilcox.test(ilisi_result$value, ilisi_result$variable, p.adj = "BH")
    table_pvalue <- test_wilcoxon$p.value
    table_pvalue_round <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'<0.001',round(table_pvalue,3)))
    table_pvalue_round <- data.frame(table_pvalue_round, check.names = FALSE) 
    write.xlsx(x = table_pvalue_round, file = paste0(save_dir, dataset2,'/statistical_test/iLISI_wilcoxon_test.xlsx'),showNA = F,sheetName = 'pvalue')
    table_pvalue2 <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'***',ifelse(table_pvalue<0.01,'**',ifelse(table_pvalue<0.05,'*','NS'))))
    table_pvalue2 <- data.frame(table_pvalue2, check.names = FALSE)
    write.xlsx(x = table_pvalue2,file = paste0(save_dir, dataset2,'/statistical_test/iLISI_wilcoxon_test.xlsx'),showNA = F,sheetName = 'pvalue2', append=TRUE)
  
    ###################
    # cLISI
    if(dataset!='dataset9'){
      print('cLISI')
      clisi_result <- read.csv(paste0(save_dir, dataset2,'/cLISI_summary.csv'), check.names = F)
      clisi_result <- melt(clisi_result)
      test_wilcoxon <- pairwise.wilcox.test(clisi_result$value, clisi_result$variable, p.adj = "BH")
      table_pvalue <- test_wilcoxon$p.value
      table_pvalue_round <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'<0.001',round(table_pvalue,3)))
      table_pvalue_round <- data.frame(table_pvalue_round, check.names = FALSE) 
      
      write.xlsx(x = table_pvalue_round, file = paste0(save_dir, dataset2,'/statistical_test/cLISI_wilcoxon_test.xlsx'),showNA = F,sheetName = 'pvalue')
      table_pvalue2 <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'***',ifelse(table_pvalue<0.01,'**',ifelse(table_pvalue<0.05,'*','NS'))))
      table_pvalue2 <- data.frame(table_pvalue2, check.names = FALSE)
      write.xlsx(x = table_pvalue2,file = paste0(save_dir, dataset2,'/statistical_test/cLISI_wilcoxon_test.xlsx'),showNA = F,sheetName = 'pvalue2', append=TRUE)
    }
  }
  # ###################
  # # ARI TEST
  # ###################
  if(ari_stat_test){
    print('ARI')

    ari_result <- read.table(paste0('manuscript_results/nicole_eval/ARI_',dataset,'/ARISampled_CT/AllValues_ARISampled_CTARI_',dataset,'.txt'), head=T, check.names = F)

    ari_result$use_case <- relevel(ari_result$use_case, ref = "LIGER")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "scMerge")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "ZINB-WaVE")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "MMD-ResNet")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "Scanorama")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "scGen")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "limma")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "ComBat")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "MNN_Correct")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "fastMNN")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "Harmony")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "Seurat_3")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "Seurat_2")
    ari_result$use_case <- relevel(ari_result$use_case, ref = "Raw")

    if(dataset=='dataset6'){
      ari_result$ari_batch <- (ari_result$ari_batchfirst + ari_result$ari_batchsecond)/2
      ari_result_batch <- ari_result
      # ari_result_batch <- data.frame(use_case=rep(ari_result$use_case,2), ari_batch=c(ari_result$ari_batchfirst,ari_result$ari_batchsecond))
    } else {
      ari_result_batch <- ari_result
    }

    test_wilcoxon_batch <- pairwise.wilcox.test(ari_result_batch$ari_batch, ari_result_batch$use_case, p.adj = "BH")
    table_pvalue <- test_wilcoxon_batch$p.value
    table_pvalue_round <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'<0.001',round(table_pvalue,3)))
    table_pvalue_round <- data.frame(table_pvalue_round, check.names = FALSE) 
    write.xlsx(x = table_pvalue_round, file = paste0(save_dir, dataset2, '/statistical_test/ARI_batch_wilcoxon_test_v2.xlsx'),showNA = F,sheetName = 'pvalue')
    table_pvalue2 <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'***',ifelse(table_pvalue<0.01,'**',ifelse(table_pvalue<0.05,'*','NS'))))
    table_pvalue2 <- data.frame(table_pvalue2, check.names = FALSE)
    write.xlsx(x = table_pvalue2,file = paste0(save_dir, dataset2, '/statistical_test/ARI_batch_wilcoxon_test_v2.xlsx'),showNA = F,sheetName = 'pvalue2', append=TRUE)

    if(dataset!='dataset9'){
      test_wilcoxon_celltype <- pairwise.wilcox.test(ari_result$ari_celltype, ari_result$use_case, p.adj = "BH")
      table_pvalue <- test_wilcoxon_celltype$p.value
      table_pvalue_round <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'<0.001',round(table_pvalue,3)))
      table_pvalue_round <- data.frame(table_pvalue_round, check.names = FALSE) 
      write.xlsx(x = table_pvalue_round, file = paste0(save_dir, dataset2,'/statistical_test/ARI_celltype_wilcoxon_test_v2.xlsx'),showNA = F,sheetName = 'pvalue')
      table_pvalue2 <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'***',ifelse(table_pvalue<0.01,'**',ifelse(table_pvalue<0.05,'*','NS'))))
      table_pvalue2 <- data.frame(table_pvalue2, check.names = FALSE)
      write.xlsx(x = table_pvalue2,file = paste0(save_dir ,dataset2, '/statistical_test/ARI_celltype_wilcoxon_test_v2.xlsx'),showNA = F,sheetName = 'pvalue2', append=TRUE)
    }
  }
  ###################
  # ASW
  ###################
  if(asw_stat_test){
    print('ASW')
    
    if(dataset=='dataset9'){
      # if(file.exists(paste0('xiaomeng/ASW/',dataset,'_.*\\ASW_batch.txt'))){
      lsfn <- list.files('xiaomeng/ASW',pattern = paste0(dataset,'_.*\\ASW_batch.txt'),full.names = TRUE)
      if(length(lsfn)!=0){
        asw_result_list <- lapply(lsfn, function(i){
            data <- read.table(i,head=T,row.names=1, check.names = F)
            colnames(data) <- c('batch','row_name')
            return(data)
        })
      }
    } else {
      lsfn <- list.files('xiaomeng/ASW',pattern = paste0(dataset,'_.*\\ASW.txt'),full.names = TRUE)
      
      if(length(lsfn)!=0){
        asw_result_list <- lapply(lsfn, function(i){
          if(!grepl('scgen_without_cell_information',i)){
            data <- read.table(i, head=T, row.names=1, check.names = F)
            colnames(data) <- c('batch','celltype','row_name')
            return(data)
          }
        })
      }  
    }

    asw_result <- do.call(rbind, asw_result_list)
    asw_result$row_name2[1:2]
    asw_result$row_name2 <- ifelse(asw_result$row_name=="classicMNN",'MNN_Correct',
                            ifelse(asw_result$row_name=="seurat3",'Seurat_3',
                                   ifelse(asw_result$row_name=="seurat2",'Seurat_2',
                                          ifelse(asw_result$row_name=="combat",'ComBat',
                                                 ifelse(asw_result$row_name=="liger",'LIGER',
                                                        ifelse(asw_result$row_name=="scmerge",'scMerge',
                                                               ifelse(asw_result$row_name=="zinbwave",'ZINB-WaVE',
                                                                      ifelse(asw_result$row_name=="resnet",'MMD-ResNet',
                                                                             ifelse(asw_result$row_name=="scanorama",'Scanorama',
                                                                                    ifelse(asw_result$row_name=="harmony",'Harmony',
                                                                                           ifelse(asw_result$row_name=="scgen",'scGen',
                                                                                                  ifelse(asw_result$row_name=="limma",'limma',
                                                                                                         ifelse(asw_result$row_name=="fastMNN",'fastMNN',
                                                                                                                ifelse(asw_result$row_name=="raw",'Raw',NA))))))))))))))
     

    
    asw_result$row_name2 <- factor(asw_result$row_name2 , ordered = FALSE)
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "LIGER")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "scMerge")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "ZINB-WaVE")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "MMD-ResNet")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "Scanorama")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "scGen")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "limma")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "ComBat")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "MNN_Correct")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "fastMNN")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "Harmony")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "Seurat_3")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "Seurat_2")
    asw_result$row_name2 <- relevel(asw_result$row_name2, ref = "Raw")
    
    test_wilcoxon <- pairwise.wilcox.test(asw_result$batch, asw_result$row_name2, p.adj = "BH")
    
    table_pvalue <- test_wilcoxon$p.value
    table_pvalue_round <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'<0.001',round(table_pvalue,3)))
    table_pvalue_round <- data.frame(table_pvalue_round, check.names = FALSE) 
    write.xlsx(x = table_pvalue_round, file = paste0(save_dir, dataset2,'/statistical_test/ASW_batch_wilcoxon_test.xlsx'),showNA = F,sheetName = 'pvalue')
    table_pvalue2 <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'***',ifelse(table_pvalue<0.01,'**',ifelse(table_pvalue<0.05,'*','NS'))))
    table_pvalue2 <- data.frame(table_pvalue2, check.names = FALSE)
    write.xlsx(x = table_pvalue2,file = paste0(save_dir, dataset2,'/statistical_test/ASW_batch_wilcoxon_test.xlsx'),showNA = F,sheetName = 'pvalue2', append=TRUE)

    if(dataset!='dataset9'){ # dataset 9 contain batch info only
      test_wilcoxon <- pairwise.wilcox.test(asw_result$celltype, asw_result$row_name2, p.adj = "BH")
      table_pvalue <- test_wilcoxon$p.value
      table_pvalue_round <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'<0.001',round(table_pvalue,3)))
      table_pvalue_round <- data.frame(table_pvalue_round, check.names = FALSE) 
      write.xlsx(x = table_pvalue_round, file = paste0(save_dir, dataset2,'/statistical_test/ASW_celltype_wilcoxon_test.xlsx'),showNA = F,sheetName = 'pvalue')
      table_pvalue2 <- ifelse(is.na(table_pvalue),NA, ifelse(table_pvalue<0.001,'***',ifelse(table_pvalue<0.01,'**',ifelse(table_pvalue<0.05,'*','NS'))))
      table_pvalue2 <- data.frame(table_pvalue2, check.names = FALSE)
      write.xlsx(x = table_pvalue2,file = paste0(save_dir, dataset2,'/statistical_test/ASW_celltype_wilcoxon_test.xlsx'),showNA = F,sheetName = 'pvalue2', append=TRUE)
    }
  }  
})

