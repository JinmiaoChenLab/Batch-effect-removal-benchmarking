##
correct_name <- function(input){
  colnames(input)[grep("cell_type",colnames(input))] = "celltype"
  colnames(input)[grep("CellType",colnames(input))] = "celltype"
  
  colnames(input)[grep("UMAP1",colnames(input))] = "umap1"
  colnames(input)[grep("UMAP2",colnames(input))] = "umap2"
  
  if(1 %in% input$batch){
    input$batch = paste0("Batch",input$batch)
  }
  return(input)
}


plot_function <- function(plot_df, xstring, ystring, plottype, plottitle, colorcode, xlabel='umap 1', ylabel='umap 2') {
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
  ) + scale_fill_manual(values = colorcode) + theme(axis.title.x = element_text( size = 20),axis.title.y = element_text( size = 20))
  # + scale_color_brewer(palette=marker_palette)
  
  return(plotobj)
}


create_plot_tsne <- function(input=list(raw_df=raw_df,zinbwave_df=zinbwave_df,scmerge_df=scmerge_df,seurat2_df=seurat2_df,
                                        harmony_df=harmony_df,scanorama_df=scanorama_df,scgen_df=scgen_df,resnet_df=resnet_df,
                                        fastMNN_df=fastMNN_df,seurat3_df=seurat3_df,correctMNN_df=correctMNN_df,
                                        combat_df=combat_df,limma_df=limma_df,liger_df=liger_df)){
  raw_df=input$raw_df
  zinbwave_df=input$zinbwave_df
  scmerge_df=input$scmerge_df
  seurat2_df=input$seurat2_df
  harmony_df=input$harmony_df
  scanorama_df=input$scanorama_df
  scgen_df=input$scgen_df
  resnet_df=input$resnet_df
  fastMNN_df=input$fastMNN_df
  seurat3_df=input$seurat3_df
  correctMNN_df=input$correctMNN_df
  combat_df=input$combat_df
  limma_df=input$limma_df
  liger_df=input$liger_df
  
  # zinbwave
  zinbwave_df$batch <- as.factor(zinbwave_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  zinbwave_df$celltype <- as.factor(zinbwave_df$celltype)
  plot_df = zinbwave_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'ZINB-WaVE'
  zinbwave_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  zinbwave_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = 'tSNE 1', ylabel = '')
  
  
  # scmerge
  scmerge_df$batch <- as.factor(scmerge_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  scmerge_df$celltype <- as.factor(scmerge_df$celltype)
  plot_df = scmerge_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'scMerge'
  scmerge_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  scmerge_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '', ylabel = '')
  
  
  # Liger
  liger_df$batch <- as.factor(liger_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  liger_df$celltype <- as.factor(liger_df$celltype)
  plot_df = liger_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'LIGER'
  liger_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  liger_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '', ylabel = '')
  
  
  # Combat
  combat_df$batch <- as.factor(combat_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  combat_df$celltype <- as.factor(combat_df$celltype)
  plot_df = combat_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'Combat'
  combat_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  combat_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '', ylabel = '')
  
  
  # limma
  limma_df$batch <- as.factor(limma_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  limma_df$celltype <- as.factor(limma_df$celltype)
  plot_df = limma_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'Limma'
  limma_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  limma_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '', ylabel = '')
  
  
  # Seurat2
  seurat2_df$batch <- as.factor(seurat2_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  seurat2_df$celltype <- as.factor(seurat2_df$celltype)
  plot_df = seurat2_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'Seurat2'
  seurat2_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  seurat2_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')
  
  
  # Harmony
  harmony_df$batch <- as.factor(harmony_df$batch)
  harmony_df$cell_type <- as.factor(harmony_df$celltype)
  plot_df = harmony_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'Harmony'
  harmony_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel = '',ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  harmony_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '',ylabel = '')
  
  
  # Scanorama
  scanorama_df$batch <- as.factor(scanorama_df$batch)
  scanorama_df$celltype <- as.factor(scanorama_df$celltype)
  plot_df = scanorama_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'Scanorama'
  #plottitle = 'Scanorama Batch'
  scanorama_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  #plottitle = 'Scanorama Cell Type'
  scanorama_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')
  
  
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
  scgen_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = 'tSNE 2')
  plottype = 'celltype'
  #plottitle = 'scGen Cell Type'
  plottitle = ''
  scgen_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')
  
  
  # ResNet
  resnet_df$batch <- as.factor(resnet_df$batch)
  resnet_df$celltype <- as.factor(resnet_df$celltype)
  plot_df = resnet_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'ResNet'
  #plottitle = 'ResNet Batch'
  resnet_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel='', ylabel='')
  plottype = 'celltype'
  #plottitle = 'ResNet Cell Type'
  plottitle = ''
  resnet_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel='', ylabel='')
  
  
  # fastMNN
  fastMNN_df$batch <- as.factor(fastMNN_df$batch)
  fastMNN_df$celltype <- as.factor(fastMNN_df$celltype)
  plot_df = fastMNN_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'fastMNN'
  #plottitle = 'fastMNN Batch'
  fastmnn_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel = '',ylabel = '')
  plottype = 'celltype'
  #plottitle = 'fastMNN Cell Type'
  plottitle = ''
  fastmnn_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '',ylabel = '')
  
  
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
  seurat3_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
  plottype = 'celltype'
  #plottitle = 'Seurat3 Cell Type'
  plottitle = ''
  seurat3_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')
  
  
  # classic MNN
  correctMNN_df$batch <- as.factor(correctMNN_df$batch)
  correctMNN_df$celltype <- as.factor(correctMNN_df$celltype)
  plot_df = correctMNN_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'MNNCorrect'
  #plottitle = 'correctMNN Batch'
  correctmnn_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
  plottype = 'celltype'
  #plottitle = 'correctMNN Cell Type'
  plottitle = ''
  correctmnn_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')
  
  
  # filtered data
  raw_df$batch <- as.factor(raw_df$batch)
  raw_df$celltype <- as.factor(raw_df$celltype)
  plot_df = raw_df
  xstring = 'tSNE1'
  ystring = 'tSNE2'
  plottype = 'batch'
  plottitle = 'Raw'
  #plottitle = 'Raw Data Batch'
  raw_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel = '',ylabel = 'tSNE 2')
  plottype = 'celltype'
  #plottitle = 'Raw Data Cell Type'
  plottitle = ''
  raw_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')
  
}

create_plot_umap <- function(input=list(raw_df=raw_df,zinbwave_df=zinbwave_df,scmerge_df=scmerge_df,seurat2_df=seurat2_df,
                             harmony_df=harmony_df,scanorama_df=scanorama_df,scgen_df=scgen_df,resnet_df=resnet_df,
                             fastMNN_df=fastMNN_df,seurat3_df=seurat3_df,correctMNN_df=correctMNN_df,
                             combat_df=combat_df,limma_df=limma_df,liger_df=liger_df,bbknn_df=bbknn_df)){
  raw_df=input$raw_df
  zinbwave_df=input$zinbwave_df
  scmerge_df=input$scmerge_df
  seurat2_df=input$seurat2_df
  harmony_df=input$harmony_df
  scanorama_df=input$scanorama_df
  scgen_df=input$scgen_df
  resnet_df=input$resnet_df
  fastMNN_df=input$fastMNN_df
  seurat3_df=input$seurat3_df
  correctMNN_df=input$correctMNN_df
  combat_df=input$combat_df
  limma_df=input$limma_df
  liger_df=input$liger_df
  bbknn_df=input$bbknn_df
  
  # filtered data
  raw_df$batch <- as.factor(raw_df$batch)
  raw_df$celltype <- as.factor(raw_df$celltype)
  plot_df = raw_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'Raw'
  #plottitle = 'Raw Data Batch'
  raw_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel = '')
  plottype = 'celltype'
  #plottitle = 'Raw Data Cell Type'
  plottitle = ''
  raw_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')
  
  # zinbwave
  zinbwave_df$batch <- as.factor(zinbwave_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  zinbwave_df$celltype <- as.factor(zinbwave_df$celltype)
  plot_df = zinbwave_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'ZINB-WaVE'
  zinbwave_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  zinbwave_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, ylabel = '')
  
  
  # scmerge
  scmerge_df$batch <- as.factor(scmerge_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  scmerge_df$celltype <- as.factor(scmerge_df$celltype)
  plot_df = scmerge_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'scMerge'
  scmerge_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  scmerge_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, ylabel = '',xlabel = '')
  
  
  # Liger
  liger_df$batch <- as.factor(liger_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  liger_df$celltype <- as.factor(liger_df$celltype)
  plot_df = liger_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'LIGER'
  liger_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  liger_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, ylabel = '',xlabel = '')
  
  
  # Combat
  combat_df$batch <- as.factor(combat_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  combat_df$celltype <- as.factor(combat_df$celltype)
  plot_df = combat_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'Combat'
  combat_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  combat_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '', ylabel = '')
  
  
  # limma
  limma_df$batch <- as.factor(limma_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  limma_df$celltype <- as.factor(limma_df$celltype)
  plot_df = limma_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'Limma'
  limma_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '', ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  limma_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '', ylabel = '')
  
  
  # Seurat2
  seurat2_df$batch <- as.factor(seurat2_df$batch) #only for seurat, to convert the numerical labels to factors (same as the output of other methods)
  seurat2_df$celltype <- as.factor(seurat2_df$celltype)
  plot_df = seurat2_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'Seurat2'
  seurat2_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  seurat2_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')
  
  
  # Harmony
  harmony_df$batch <- as.factor(harmony_df$batch)
  harmony_df$cell_type <- as.factor(harmony_df$celltype)
  plot_df = harmony_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'Harmony'
  harmony_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel = '',ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  harmony_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '',ylabel = '')
  
  
  # Scanorama
  scanorama_df$batch <- as.factor(scanorama_df$batch)
  scanorama_df$celltype <- as.factor(scanorama_df$celltype)
  plot_df = scanorama_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'Scanorama'
  #plottitle = 'Scanorama Batch'
  scanorama_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  #plottitle = 'Scanorama Cell Type'
  scanorama_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,ylabel = '',xlabel = '')
  
  
  # BBKNN
  bbknn_df$batch <- as.factor(bbknn_df$batch)
  bbknn_df$celltype <- as.factor(bbknn_df$celltype)
  plot_df = bbknn_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'BBKNN'
  #plottitle = 'BBKNN Batch'
  bbknn_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
  plottype = 'celltype'
  plottitle = ''
  #plottitle = 'BBKNN Cell Type'
  bbknn_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,ylabel = '',xlabel = '')
  
  
  # scGen
  scgen_df$batch <- as.factor(scgen_df$batch)
  scgen_df$celltype <- as.factor(scgen_df$celltype)
  plot_df = scgen_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'scGen'
  #plottitle = 'scGen Batch'
  scgen_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '')
  plottype = 'celltype'
  #plottitle = 'scGen Cell Type'
  plottitle = ''
  scgen_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '', ylabel = '')
  
  
  # ResNet
  resnet_df$batch <- as.factor(resnet_df$batch)
  resnet_df$celltype <- as.factor(resnet_df$celltype)
  plot_df = resnet_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'ResNet'
  #plottitle = 'ResNet Batch'
  resnet_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel='', ylabel='')
  plottype = 'celltype'
  #plottitle = 'ResNet Cell Type'
  plottitle = ''
  resnet_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel='', ylabel='')
  
  
  # fastMNN
  fastMNN_df$batch <- as.factor(fastMNN_df$batch)
  fastMNN_df$celltype <- as.factor(fastMNN_df$celltype)
  plot_df = fastMNN_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'fastMNN'
  #plottitle = 'fastMNN Batch'
  fastmnn_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch, xlabel = '',ylabel = '')
  plottype = 'celltype'
  #plottitle = 'fastMNN Cell Type'
  plottitle = ''
  fastmnn_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype, xlabel = '',ylabel = '')
  
  
  # Seurat 3
  seurat3_df$batch <- as.factor(seurat3_df$batch)
  seurat3_df$celltype <- as.factor(seurat3_df$celltype)
  # unique(seurat3_df$batch)
  plot_df = seurat3_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'Seurat3'
  #plottitle = 'Seurat3 Batch'
  seurat3_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
  plottype = 'celltype'
  #plottitle = 'Seurat3 Cell Type'
  plottitle = ''
  seurat3_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')
  
  
  # classic MNN
  correctMNN_df$batch <- as.factor(correctMNN_df$batch)
  correctMNN_df$celltype <- as.factor(correctMNN_df$celltype)
  plot_df = correctMNN_df
  xstring = 'umap1'
  ystring = 'umap2'
  plottype = 'batch'
  plottitle = 'MNNCorrect'
  #plottitle = 'correctMNN Batch'
  correctmnn_batch <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_batch,xlabel = '',ylabel = '')
  plottype = 'celltype'
  #plottitle = 'correctMNN Cell Type'
  plottitle = ''
  correctmnn_cellt <<- plot_function(plot_df, xstring, ystring, plottype, plottitle, colorcode_celltype,xlabel = '',ylabel = '')

}




