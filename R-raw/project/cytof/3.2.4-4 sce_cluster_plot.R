#!/opt/conda/bin/Rscript

#' SingleCellExperiment对象操作之：分群后绘图
#' 
#' @param sceObject sce对象的变量名
#' @param saveplotpath 保存绘图的路径
#' @param features 分群的抗体，features = "type" 表示谱系标签/ "state" 表示功能标签/ NULL 无指定,默认NULL
#' @param by_row 用来进行行聚类的数据，sample_id,cluster_id，其中的一个默认cluster_id
#' @param maxK 对应分群时使用的maxK，即分群的最大群数，默认20
#' @param hm1 分群的抗体，绘图时等同于features，可选"type","state",NULL
#' @param hm2 绘制复合表达热图时可以指定单个抗体，例如：CD4
#' @param width 图片宽度
#' @param height 图片长度
#' @param dpi 图片清晰度
#' @param ...
#' 
#' @export
sce_cluster_plot <- function(sceObject,
                             saveplotpath = "3.Cluster/",
                             features = NULL,
                             by_row = "cluster",
                             maxK = 20,
                             hm1 = NULL,
                             hm2 = "CD4",
                             bar_fontsize = 8,
                             color = SelectColors("corona"),
                             samplecolor = SelectColors("customecol2"),
                             ...){

  suppressMessages(library("CATALYST"))
  suppressMessages(library("cytofWorkflow"))
  suppressMessages(library("ggplot2"))
  suppressMessages(library("ComplexHeatmap"))
  suppressMessages(library("cowplot"))
  suppressMessages(library("reshape2"))
  sceObject <- readdata(filename = sceObject)
  
  # 密度图 25*17
  picture2 <- plotClusterExprs(sceObject, k = paste0("meta", maxK), features = NULL) + ylab("cluster")
  ggplotsave(plot = picture2, 
             savepath = saveplotpath, 
             mapname = "sce_cluster_ClusterExprs", 
             width = 25, height = 10,...)
  
  ns <- table(cluster_id = cluster_ids(sceObject, paste0('meta', maxK)), 
              sample_id = sample_ids(sceObject))
  fq <- prop.table(ns, 2) * 100
  df3 <- as.data.frame(fq)
  df3.1 <- dcast(df3, cluster_id~sample_id)
  savexlsx1(data = df3.1,
            filename = paste0(saveplotpath, "/sce_cluster_ratio_info.xlsx"),
            sheet = "sce_cluster_ratio_info")
  
  # 堆叠柱状图
  plot_arg <- list(small_sample = scale_x_discrete(expand = c(0, 0.5)), 
                   middle_sample = scale_x_discrete(expand = c(0, 1)),
                   big_sample = scale_x_discrete(expand = c(0, 2)))
  if(length(table(sceObject@colData$sample_id)) < 3){
    plot_arg_needed <- plot_arg$small_sample
  }
  if(length(table(sceObject@colData$sample_id)) > 2 & length(table(sceObject@colData$sample_id)) < 6){
    plot_arg_needed <- plot_arg$middle_sample
  }
  if(length(table(sceObject@colData$sample_id)) > 5 & length(table(sceObject@colData$sample_id)) < 9){
    plot_arg_needed <- plot_arg$big_sample
  }
  if(length(table(sceObject@colData$sample_id)) > 8){
    plot_arg_needed <- NULL
  }
  picture3 <- plotAbundances(sceObject, k = paste0("meta", maxK), by = "sample") + 
    plot_arg_needed+
    scale_fill_manual(values = samplecolor)+ 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
          axis.text = element_text(colour = "black", size = bar_fontsize))
  
  ggplotsave(plot = picture3,
             savepath = saveplotpath, mapname = "sce_cluster_ratio", 
             width = 8, height = 8,)
  
  # 箱线图 16*9
  picture4 <- plotAbundances(sceObject, k = paste0("meta", maxK), by = "cluster") + 
    scale_color_manual(values = color) + 
    scale_fill_manual(values = c(rep("#FFFFFF", length(table(sceObject$group))))) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) + labs(color = "group")
  
  ggplotsave(plot = picture4, 
             savepath = saveplotpath, 
             mapname = "sce_cluster_box", 
             width = 16, height = 7,...)
  #write.xlsx(data.frame(picture4$data), paste0(saveplotpath, "box_info.xlsx"), rowNames = F, colNames = T)
  
  plotfile(savepath = saveplotpath,
           mapname = "sce_cluster_ExprHeatmap",
           width = 10,height = 8,...)
  
  picture5 <- plotExprHeatmap(sceObject, features = features, # features = "type"/ "state"/ NULL
                              by = by_row, k = paste0("meta", maxK), # 根据选定的maxK来调整 不能超过maxK
                              row_anno = "cluster",scale = "last",
                              bars = TRUE, perc = TRUE)
  picture5@left_annotation@anno_list[["cluster_id"]]@fun@var_env[["color_mapping"]]@full_col[1:length(table(sceObject$cluster))] <- samplecolor[1:length(table(sceObject$cluster))]
  picture5@left_annotation@anno_list[["cluster_id"]]@color_mapping@full_col[1:length(table(sceObject$cluster))] <- samplecolor[1:length(table(sceObject$cluster))]
  draw(picture5)
  plotsave()

  plotfile(savepath = saveplotpath,
           mapname = "sce_cluster_ExprHeatmap_no_clust",
           width = 10,height = 8,...)
  
  picture6 <- plotExprHeatmap(sceObject, features = features, # features = "type"/ "state"/ NULL
                              by = by_row, k = paste0("meta", maxK), # 根据选定的maxK来调整 不能超过maxK
                              row_anno = "cluster",scale = "last", row_clust = F, col_clust = F,
                              bars = TRUE, perc = TRUE)
  picture6@left_annotation@anno_list[["cluster_id"]]@fun@var_env[["color_mapping"]]@full_col[1:length(table(sceObject$cluster))] <- samplecolor[1:length(table(sceObject$cluster))]
  picture6@left_annotation@anno_list[["cluster_id"]]@color_mapping@full_col[1:length(table(sceObject$cluster))] <- samplecolor[1:length(table(sceObject$cluster))]
  draw(picture6)
  plotsave()

  df5 <- as.data.frame(picture5@matrix)
  rownames(df5) <- paste0("cluster", rownames(df5))
  
  savexlsx3(data = df5,
            filename = paste0(saveplotpath, "/heatmap_info.xlsx"),
            sheet = "heatmap_info")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){

  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-sce","--sceObject",default = "rawdata/flowSet_object_sce_cluster.rds", help = "读取sce对象的rds文件名")
  parser$add_argument("-by","--by_row",default = "cluster", help = "用来进行行聚类的数据，sample_id,cluster_id，其中的一个;默认cluster_id")
  parser$add_argument("-max","--maxK",default = 20, type = "integer", help = "对应分群时使用的maxK，即分群的最大群数，默认20")
  parser$add_argument("-h1","--hm1",default = NULL, help = "分群的抗体，绘图时等同于features，可选'type','state',NULL")
  parser$add_argument("-h2","--hm2",default = "CD4", help = "绘制复合表达热图时可以指定单个抗体，例如：CD4")
  parser$add_argument("-fe","--features",default = NULL, help = "抗体类型，features = 'type' 表示谱系标签/ 'state' 表示状态标签/ NULL 无指定,默认NULL")
  
  parser$add_argument("-sp","--saveplotpath",default = "3.Cluster/", help = "保存图片的目录")
  parser$add_argument("-i","--dpi", default = 300, type = "integer", help = "图片清晰度")
  parser$add_argument("-fa","--family", default = "sans", help = "图片字体")
  parser$add_argument("-font","--bar_fontsize", default = 8, type = "integer", help = "图片字体大小")
  
  args <- parser$parse_args()
  
  clustplotargs <- do.call(what = sce_cluster_plot, args = args)
}
## 注意：每个目录和文件必须严格存在！rawdata是固定的下机数据存放目录 flowSet_object.rds是固定的文件名
