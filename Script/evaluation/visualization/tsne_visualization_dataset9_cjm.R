rm(list=ls())

source("functions.R")

#libraries
library(ggplot2)
library(gridExtra)

# Script for visualization 
this.dir <- "T:/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/manuscript_results/visualization/"
setwd(this.dir)

input_path = "T:/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/dataset9/"
raw_path = "T:/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/dataset9/"


dataset = "dataset9"
base_name = 'dataset9'

###############################
#load data

seurat2_df <- read.csv(paste0(input_path,dataset,"_seurat2_tsne.csv"), head=T, row.names = 1, check.names = FALSE)
seurat2_df=correct_name(seurat2_df)

harmony_df <- read.csv(paste0(input_path,dataset,"_harmony_tsne.csv"), head=T, row.names = 1, check.names = FALSE)
harmony_df=correct_name(harmony_df)
  
scanorama_df <- read.csv(paste0(input_path,dataset,"_scanorama_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
scanorama_df=correct_name(scanorama_df)
# bbknn_df <- read.csv("bbknn_tsne.csv", row.names = 1, check.names = FALSE)

scgen_df <- read.csv(paste0(input_path,dataset,"_scgen_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
scgen_df=correct_name(scgen_df)

resnet_df <- read.csv(paste0(input_path,dataset,"_resnet_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
resnet_df=correct_name(resnet_df)

#raw_df <- read.csv(paste0(raw_path,"rawdata_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
#raw_df <- read.table(paste0(raw_path,"filtered_tsne.txt"),head=T, row.names = 1, check.names = FALSE)
raw_df <- read.csv(paste0(raw_path,dataset,"_raw_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
raw_df=correct_name(raw_df)
#colnames(raw_df)[grep("tSNE_1",colnames(raw_df))] = "tSNE1"
#colnames(raw_df)[grep("tSNE_2",colnames(raw_df))] = "tSNE2"


fastMNN_df <- read.csv(paste0(input_path,dataset,"_fastMNN_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
fastMNN_df=correct_name(fastMNN_df)

seurat3_df <- read.csv(paste0(input_path,dataset,"_seurat3_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
seurat3_df=correct_name(seurat3_df)

correctMNN_df <- read.csv(paste0(input_path,dataset,"_classicMNN_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
correctMNN_df=correct_name(correctMNN_df)

combat_df <- read.csv(paste0(input_path,dataset,"_combat_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
combat_df=correct_name(combat_df)

limma_df <- read.csv(paste0(input_path,dataset,"_limma_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
limma_df=correct_name(limma_df)

liger_df <- read.csv(paste0(input_path,dataset,"_liger_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
liger_df=correct_name(liger_df)

scmerge_df <- read.csv(paste0(input_path,dataset,"_scmerge_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
scmerge_df=correct_name(scmerge_df)

zinbwave_df <- read.csv(paste0(input_path,dataset,"_zinbwave_tsne.csv"),head=T, row.names = 1, check.names = FALSE)
zinbwave_df=correct_name(zinbwave_df)
# seurat2_df <- read.table("./seurat2_multicca/seuratCCA_tsne.txt", head=T, row.names = 1, check.names = FALSE)
# harmony_df <- read.table("./harmony/harmony_tsne.txt", head=T, row.names = 1, check.names = FALSE)
# scanorama_df <- read.csv("./scanorama/scanorama_tsne.csv", row.names = 1, check.names = FALSE)
# # bbknn_df <- read.csv("./bbknn/bbknn_tsne.csv", row.names = 1, check.names = FALSE)
# scgen_df <- read.csv("./scGen/scGen_tsne.csv", row.names = 1, check.names = FALSE)
# resnet_df <- read.csv("./resnet/resnet_tsne.csv", row.names = 1, check.names = FALSE)
# raw_df <- read.table("./raw/filtered_tsne.txt", head=T, row.names = 1, check.names = FALSE)
# fastMNN_df <- read.table("./fastMNN/fastMNN_tsne.txt", head=T, row.names = 1, check.names = FALSE)
# seurat3_df <- read.table("./seurat3_integration/seurat3_tsne.txt", head=T, row.names = 1, check.names = FALSE)
# correctMNN_df <- read.csv("./MNNCorrect/mnncorrect_scanpy_tsne_cellatlas.csv", row.names = 1, check.names = FALSE)
# combat_df <- read.table("./combat/combat_corrected_dataset2_tsne.txt", head=T, row.names = 1, check.names = FALSE)
# # limma_df <- read.table("./limma/limma_tsne.txt", head=T, row.names = 1, check.names = FALSE)
# liger_df <- read.table("./liger/liger_tsne.txt", head=T, row.names = 1, check.names = FALSE)
# scmerge_df <- read.csv("./scMerge/scmerge_tsne.csv", row.names = 1, check.names = FALSE)
# # zinbwave_df <- read.table("./zinbwave/zinbwave_tsne.txt", head=T, row.names = 1, check.names = FALSE)

#for checking purposes
#View(head(scanorama_df))  
#summary(as.factor(seurat2_df$batch))

###############################
#define plot funs

#colour palettes used: 
#http://www.sthda.com/english/wiki/colors-in-r



# element_line(color = "black", size = 0.5, linetype = "solid")
# +
#   scale_y_continuous( breaks=c(ymin,ymax) ) + 
#   scale_x_continuous( breaks=c(xmin,xmax) ) + 
#   coord_cartesian(xlim =c(xmin,xmax), ylim = c(ymin,ymax))
# marker_palette
plot_function <- function(plot_df, xstring, ystring, plottype, plottitle, colorcode, xlabel='tSNE 1', ylabel='tSNE 2') {
  plotobj <- ggplot(plot_df, aes_string(x = xstring, y = ystring, colour = plottype)) + 
    geom_point(shape=19, size=0.7) + theme_bw() #alpha = 0.6
  #ggplot_build(plotobj)$data
  
  xmax = max(plot_df[[xstring]])
  xmin = min(plot_df[[xstring]])
  ymax = max(plot_df[[ystring]])
  ymin = min(plot_df[[ystring]])
  xmax = ceiling(xmax)
  xmin = floor(xmin)
  ymax = ceiling(ymax)
  ymin = floor(ymin)
  
  if(plottitle==''){
    plotobj <- plotobj + labs(x=xlabel,y=ylabel) 
  } else {
    plotobj <- plotobj + labs(x=xlabel,y=ylabel,title=plottitle) 
  }
  
  plotobj <- plotobj + theme(legend.title = element_blank(), 
                             plot.title = element_text(color="black", size=20,hjust = 0.5),
                             legend.position = "none", 
                             axis.line = element_blank(), 
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                             panel.border = element_blank(),
                             axis.text = element_blank(),
                             #axis.title = element_blank(),
                             axis.ticks = element_blank()
  ) + scale_fill_manual(values = colorcode)
  # + scale_color_brewer(palette=marker_palette)
  
  return(plotobj)
}


legend_function <- function(plot_df, xstring, ystring, plottype, plottitle, colorcode, legendtitle, legendlabels=NULL) {
  #theme_bw() + 
  #alpha = 0.6
  # scale_color_brewer(palette=marker_palette, labels = legendlabels) +  legendtitle
  plotobj <- ggplot(plot_df, aes_string(x = xstring, y = ystring, colour = plottype)) +
    geom_point(shape=19, size=5) +
    scale_fill_manual(values = colorcode) + labs(colour=legendtitle)
  
  plotobj <- plotobj + theme(
    legend.title = element_text(size=20), 
    legend.text = element_text(size=15), 
    plot.title = element_text(color="black", size=20, hjust = 0.5)
  )
  #legend.key.size = unit(1.1, "cm"),
  #legend.key.width = unit(0.5,"cm"),
  return(plotobj)
}

###############################
# generate the plot lists(?)

batch_palette = "Set1"
celltype_palette = "Dark2"
library(RColorBrewer)
colourCountB = length(unique(seurat2_df$batch))
colourCountCT = length(unique(seurat2_df$celltype))
# getPaletteB = colorRampPalette(brewer.pal(8, "Set1"))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
colorcode_celltype <- getPalette(colourCountCT)
# colorcode_batch <- getPaletteB(colourCountB)
colorcode_batch <- c("#E41A1C","#268AB5")


# zinbwave
zinbwave_df$batch <- as.factor(zinbwave_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
zinbwave_df$celltype <- as.factor(zinbwave_df$celltype)
plot_df = zinbwave_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'ZINB-WaVE'
zinbwave_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
plottype = 'celltype'
plottitle = ''
zinbwave_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, ylabel = '')


# scmerge
scmerge_df$batch <- as.factor(scmerge_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
scmerge_df$celltype <- as.factor(scmerge_df$celltype)
plot_df = scmerge_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'scMerge'
scmerge_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
plottype = 'celltype'
plottitle = ''
scmerge_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, ylabel = '')


# Liger
liger_df$batch <- as.factor(liger_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
liger_df$celltype <- as.factor(liger_df$celltype)
plot_df = liger_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'LIGER'
liger_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
plottype = 'celltype'
plottitle = ''
liger_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, ylabel = '')


# Combat
combat_df$batch <- as.factor(combat_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
combat_df$celltype <- as.factor(combat_df$celltype)
plot_df = combat_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'Combat'
combat_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
plottype = 'celltype'
plottitle = ''
combat_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '', ylabel = '')


# limma
limma_df$batch <- as.factor(limma_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
limma_df$celltype <- as.factor(limma_df$celltype)
plot_df = limma_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'Limma'
limma_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
plottype = 'celltype'
plottitle = ''
limma_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '', ylabel = '')


# Seurat2
seurat2_df$batch <- as.factor(seurat2_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
seurat2_df$celltype <- as.factor(seurat2_df$celltype)
plot_df = seurat2_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'Seurat2'
seurat2_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
plottype = 'celltype'
plottitle = ''
seurat2_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')


# Harmony
harmony_df$batch <- as.factor(harmony_df$batch)
harmony_df$cell_type <- as.factor(harmony_df$celltype)
plot_df = harmony_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'Harmony'
harmony_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel = '',ylabel = '')
plottype = 'celltype'
plottitle = ''
harmony_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '',ylabel = '')


# Scanorama
scanorama_df$batch <- as.factor(scanorama_df$batch)
scanorama_df$celltype <- as.factor(scanorama_df$celltype)
plot_df = scanorama_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'Scanorama'
#plottitle = 'Scanorama Batch'
scanorama_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
plottype = 'celltype'
plottitle = ''
#plottitle = 'Scanorama Cell Type'
scanorama_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,ylabel = '')


# BBKNN
# bbknn_df$batch <- as.factor(bbknn_df$batch)
# bbknn_df$celltype <- as.factor(bbknn_df$celltype)
# plot_df = bbknn_df
# xstring = 'tSNE1'
# ystring = 'tSNE2'
# plottype = 'batch'
# plottitle = 'BBKNN'
# #plottitle = 'BBKNN Batch'
# bbknn_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
# plottype = 'celltype'
# plottitle = ''
# #plottitle = 'BBKNN Cell Type'
# bbknn_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,ylabel = '')


# scGen
scgen_df$batch <- as.factor(scgen_df$batch)
scgen_df$celltype <- as.factor(scgen_df$celltype)
plot_df = scgen_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'scGen'
#plottitle = 'scGen Batch'
scgen_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '')
plottype = 'celltype'
#plottitle = 'scGen Cell Type'
plottitle = ''
scgen_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype)


# ResNet
resnet_df$batch <- as.factor(resnet_df$batch)
resnet_df$celltype <- as.factor(resnet_df$celltype)
plot_df = resnet_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'ResNet'
#plottitle = 'ResNet Batch'
resnet_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel='', ylabel='')
plottype = 'celltype'
#plottitle = 'ResNet Cell Type'
plottitle = ''
resnet_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, ylabel='')


# fastMNN
fastMNN_df$batch <- as.factor(fastMNN_df$batch)
fastMNN_df$celltype <- as.factor(fastMNN_df$celltype)
plot_df = fastMNN_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'fastMNN'
#plottitle = 'fastMNN Batch'
fastmnn_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel = '',ylabel = '')
plottype = 'celltype'
#plottitle = 'fastMNN Cell Type'
plottitle = ''
fastmnn_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '',ylabel = '')


# Seurat 3
seurat3_df$batch <- as.factor(seurat3_df$batch)
seurat3_df$celltype <- as.factor(seurat3_df$celltype)
# unique(seurat3_df$batch)
plot_df = seurat3_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'Seurat3'
#plottitle = 'Seurat3 Batch'
seurat3_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
plottype = 'celltype'
#plottitle = 'Seurat3 Cell Type'
plottitle = ''
seurat3_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')


# classic MNN
correctMNN_df$batch <- as.factor(correctMNN_df$batch)
correctMNN_df$celltype <- as.factor(correctMNN_df$celltype)
plot_df = correctMNN_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'MNNCorrect'
#plottitle = 'correctMNN Batch'
correctmnn_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
plottype = 'celltype'
#plottitle = 'correctMNN Cell Type'
plottitle = ''
correctmnn_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')


# filtered data
raw_df$batch <- as.factor(raw_df$batch)
raw_df$celltype <- as.factor(raw_df$celltype)
plot_df = raw_df
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
plottitle = 'Raw'
#plottitle = 'Raw Data Batch'
raw_batch = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel = '')
plottype = 'celltype'
#plottitle = 'Raw Data Cell Type'
plottitle = ''
raw_cellt = plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '')

source("functions.R")
create_plot_tsne()

###############################
# plot setup

plot_df = harmony_df
plot_df$batchlabel <- ifelse(plot_df$batch==1,"Batch 1",ifelse(plot_df$batch==2,"Batch 2","Others"))
# harmony_df$batch <- as.factor(harmony_df$batch)

###
plot_df$celltype = as.character(plot_df$celltype)
plot_df$celltype[plot_df$celltype=="Plasmacytoid dendritic cell"] = "Plasmacytoid\ndendritic cell"
plot_df$celltype[plot_df$celltype=="Hematopoietic stem cell"] = "Hematopoietic\nstem cell"
plot_df$celltype[plot_df$celltype=="Monocyte_CD14"] = "Monocyte CD14"
plot_df$celltype[plot_df$celltype=="Monocyte_FCGR3A"] = "Monocyte FCGR3A"

plot_df$batch = paste0(plot_df$batch,"(",plot_df$tissue,")")
###
unique(harmony_df$batch)
t <- unique(harmony_df$cell_type)  # for verification
plottitle = ""
xstring = 'tSNE1'
ystring = 'tSNE2'
plottype = 'batch'
legendtitle = "Batch"
# legendlabels = c('Batch 1', 'Batch 2')
# marker_palette = batch_palette
p2a = legend_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, legendtitle)

plottype = 'celltype'
legendtitle = "Cell Type"
# legendlabels = c('CD141', 'CD1C', 'Double Neg', 'pDC')
# marker_palette = celltype_palette
p2b = legend_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, legendtitle)
p2b

# output figure
legend1 <- cowplot::get_legend(p2a)
legend2 <- cowplot::get_legend(p2b)
legend11 = legend1
legend21 = legend2
# #8x3
# 
# layout <- rbind(c(1,2,3,4,5,21), c(6,7,8,9,10,22),
#                 c(11,12,13,14,15,23), c(16,17,18,19,20,24))
###############################
# Final plot


#layout <- rbind(c(1,2,3,4,5,6,7,8), 
#                c(9,10,11,12,13,14,NA,16),
#                c(17,18,19,20,21,22,NA,24),
#                c(25,26,27,28,29,30,NA,32))

layout <- rbind(c(1,2,3,4,5,6,7,8), 
                c(NA,NA,NA,NA,NA,NA,NA,NA),
                c(17,18,19,20,21,22,NA,24),
                c(NA,NA,NA,NA,NA,NA,NA,NA))

legend2 <- legend1

plots <- list(raw_batch, seurat2_batch, seurat3_batch, harmony_batch, fastmnn_batch, correctmnn_batch, combat_batch, limma_batch,
              raw_cellt, seurat2_cellt, seurat3_cellt, harmony_cellt, fastmnn_cellt, correctmnn_cellt, combat_cellt, limma_cellt,
              scgen_batch, scanorama_batch, resnet_batch, zinbwave_batch, scmerge_batch, liger_batch, raw_batch, legend1,
              scgen_cellt, scanorama_cellt, resnet_cellt, zinbwave_cellt, scmerge_cellt, liger_cellt, raw_cellt, legend2)

# plots <- list(raw_batch, seurat2_batch, seurat3_batch, harmony_batch, fastmnn_batch, correctmnn_batch, combat_batch, raw_batch,
#               raw_cellt, seurat2_cellt, seurat3_cellt, harmony_cellt, fastmnn_cellt, correctmnn_cellt, combat_cellt, raw_cellt,  
#               scgen_batch, scanorama_batch, resnet_batch, raw_batch, scmerge_batch, liger_batch, raw_batch, legend1,
#               scgen_cellt, scanorama_cellt, resnet_cellt, raw_cellt, scmerge_cellt, liger_cellt, raw_cellt, legend2)

# bbknn have umap plot only, no tsne

#plotheight = 2*800
#plotwidth = 2*1000
#plotres = 2*72
#plotheight = 8 * 4
#plotwidth = 3 * 4

select_grobs <- function(lay) {
  id <- unique(c(t(lay)))
  id[!is.na(id)]
}

png(paste(base_name,"_tsne.png",sep=""),height = 2*900, width=2*1700,res = 2*72)
# pdf(paste(base_name,"_tsne.pdf",sep=""), height=plotheight, width=plotwidth)
print(grid.arrange(grobs = plots[select_grobs(layout)], layout_matrix = layout,
                   bottom=" ",right=" "))
#top="TSNE Visualization"
dev.off()

#png(paste(base_name,"_tsne.png",sep=""),height = 2*900, width=2*1700,res = 2*72)
pdf(paste(base_name,"_tsne.pdf",sep=""), height=2*900/300*2, width=2*1700/300*2)
print(grid.arrange(grobs = plots[select_grobs(layout)], layout_matrix = layout,
                   bottom=" ",right=" "))
#top="TSNE Visualization"
dev.off()

