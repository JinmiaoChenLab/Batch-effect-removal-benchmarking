source('../../ZINB_WaVE/ZINB_WaVE_analysis.R')
library(fs)
library(parallel)

nThread = 8
# base_name = basename(getwd())

files = dir('../../Data/dataset3/', '.*', full.names = T)
# # #### preprocessing ####
# 
# all_rc_and_anno_list = lapply(1:length(files), function(i){
#   if(is_file(files[i])){
#     return(NULL)
#   }
#   print(files[i])
#   data_dir = files[i]
#   # geneinfo <- read.table(paste0(data_dir,"/geneinfo.txt"), head=T, row.names = 1, check.names = FALSE)
#   cellinfo <- read.table(paste0(data_dir,"/cellinfo.txt"), head=T, row.names = 1, check.names = FALSE)
#   counts <- read.table(paste0(data_dir,"/counts.txt"), head=T, row.names = 1, check.names = FALSE)
# 
#   counts_hm <- t(counts)
#   dim(counts_hm)
#   rs = rowSums(counts_hm)
#   counts_hm = counts_hm[rs > 0, , drop = F]
#   dim(counts_hm)
#   ez_head(counts_hm)
#   colnames(counts_hm) <- paste(cellinfo[colnames(counts_hm),'Group'],colnames(counts_hm),sep="_")
#   rownames(cellinfo) <- paste(cellinfo[rownames(cellinfo),'Group'],rownames(cellinfo), sep="_")
#   colnames(cellinfo)[2:3] = c("batch", "cell_type")
#   ez_head(cellinfo)
#   cellinfo <- cellinfo[colnames(counts_hm), , drop = F]
#   res = list(counts = counts_hm, cell_anno = cellinfo)
#   return(res)
# })
# names(all_rc_and_anno_list) = basename(files)
# all_rc_and_anno_list2 = all_rc_and_anno_list[!sapply(all_rc_and_anno_list, is.null)]
# saveRDS(all_rc_and_anno_list2, "all_rc_and_anno_list.RDS")
# # #######################

all_rc_and_anno_list = fast_read_table('./all_rc_and_anno_list.RDS')

# lapply(1:length(all_rc_and_anno_list), function(i){
#   print(names(all_rc_and_anno_list)[i])
#   ZINB_WaVE_analysis(all_rc_and_anno_list[[i]]$counts
#                      , all_rc_and_anno_list[[i]]$cell_anno
#                      , base_name = names(all_rc_and_anno_list)[i]
#                      , nthread = nThread)
# })

mclapply(1:length(all_rc_and_anno_list), function(i){
  print(names(all_rc_and_anno_list)[i])
  ZINB_WaVE_analysis(all_rc_and_anno_list[[i]]$counts
                     , all_rc_and_anno_list[[i]]$cell_anno
                     , base_name = names(all_rc_and_anno_list)[i]
                     , nthread = nThread)
}, mc.cores = 8)
