#!/opt/conda/bin/Rscript

#' SingleCellExperiment对象操作之：质控绘图
#' 
#' @param sceObject sce对象的变量名
#' @param color_by 绘图的颜色以该参数区分，可选：sample_id、condition、patient_id
#' @param facet_col 密度图分面的列数
#' @param group_by 确定绘图的分组，可选：sample_id、condition、patient_id
#' @param label_by 定义图中的标签显示
#' @param scale_method 中位数标准化的方法，可选first、last、never 无指定,默认last
#' @param features 抗体类型，features = "type" 表示谱系标签/ "state" 表示状态标签/ NULL 无指定,默认NULL
#' @param saveplotpath 保存图片的目录
#' @param width 图片宽度
#' @param height 图片长度
#' @param dpi 图片清晰度
#' @param family 图片字体
#' 
#' @export
sce_QC_plot <- function(sceObject,
                        color_by = "group",
                        facet_col = 5,
                        group_by = "sample",
                        label_by ="sample",
                        scale_method = "last",
                        features = NULL,
                        saveplotpath = "1.QC/",
                        bar_fontsize = 5,
                        fontsize = 8,
                        color = SelectColors("corona"),
                        samplecolor = SelectColors("customecol2"),
                        ...){
  suppressMessages(library("ggplot2"))
  suppressMessages(library("CATALYST"))
  suppressMessages(library("cytofWorkflow"))
  suppressMessages(library("ComplexHeatmap"))
  
  sceObject <- readdata(filename = sceObject)
  
  # 样本数据展示 密度图 16*9
  picture <- plotExprs(sceObject, color_by = color_by) + 
    scale_color_manual(values = color)
  picture$facet$params$ncol <- facet_col # 分面展示的列数
  ggplotsave(plot = picture, 
             savepath = saveplotpath, 
             height = 9, width = 16,
             mapname = "sce_Per_sample_smoothed_densities_plot",
             ...) # smoothed表示信号值已经归一化
  
  savexlsx1(data = sceObject@metadata$experiment_info,
            filename = paste0(saveplotpath,"/bar_info.xlsx"),
            sheet = "bar_info")

  # 样本数量条形图Barplot showing the number of cells measured for each sample 8*8
  plot_arg <- list(small_sample = scale_x_discrete(expand = c(0, 3)), 
                   middle_sample = scale_x_discrete(expand = c(0, 2)), 
                   big_sample = scale_x_discrete(expand = c(0, 1)))
  
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
  hh1<-data.frame(table(sceObject@colData$sample))
  picture1 <- plotCounts(sceObject, group_by = group_by, color_by = color_by) + 
    plot_arg_needed +
    scale_y_continuous(limits = c(0, max(n_cells(sceObject))*1.1), expand = c(0,0)) +
    scale_fill_manual(values = color)+
    geom_text(aes(label=hh1$Freq), size=bar_fontsize/2, vjust=-0.5) + 
    theme_bw()+ 
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = bar_fontsize),
          legend.text = element_text(colour = "black", size = bar_fontsize))
    
  ggplotsave(plot = picture1,
             savepath = saveplotpath, 
             height = 6,width = 6,
             mapname = "sce_Per_sample_cellnumber_bar_plot",
             ...)
  
  # 多维尺度图 8*8
  if(length(table(sceObject@colData$sample_id)) > 2){
    picture2 <- pbMDS(sceObject, color_by = color_by, label_by = label_by) + 
      scale_color_manual(values = color)
    
    ggplotsave(plot = picture2,
               savepath = saveplotpath, 
               width = 8, height = 8,
               mapname = "sce_Per_sample_MDS_plot", ...)
    
    df2 <- data.frame(sample = picture2$data$sample, group = picture2$data$group, n_cells = picture2$data$n_cells, 
                      MDS.x = picture2$data$x, MDS.y = picture2$data$y)
    
    savexlsx1(data = df2,
              filename = paste0(saveplotpath, "/mds_info.xlsx"),
              sheet = "mds_info")
    
  }
  
  ## 非冗余丰度打分图 16*9
  # picture4 <- plotNRS(sceObject, features = features, color_by = color_by) + 
  #   theme_bw() +
  #   theme(axis.text = element_text(colour = "black", size = fontsize),
  #         legend.text = element_text(colour = "black", size = fontsize),
  #         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) + 
  #   scale_color_manual(values = color)
  # 
  # ggplotsave(plot = picture4,
  #            width = 20, height = 8,
  #            savepath = saveplotpath, 
  #            mapname = "sce_Per_sample_NRS_plot",...)
  # 
  # df4 <- picture4["data"]
  # 
  # savexlsx1(data = df4,
  #           filename = paste0(saveplotpath, "/nrs_info.xlsx"),
  #           sheet = "nrs_info")
  
  # 热图 8*8
  if(length(table(sceObject@colData$sample_id)) > 1){
    
    plotfile(savepath = saveplotpath,
             mapname = "sce_scale_ExprHeatmap.pdf",
             width = 10,height = 8,...)
    
    picture5 <- plotExprHeatmap(sceObject, scale = "last", row_anno = c("sample", "group"),
                                hm_pal = rev(hcl.colors(10, "YlGnBu")))
    picture5@left_annotation@anno_list[["sample"]]@color_mapping@full_col[1:length(table(sceObject$sample))] <- samplecolor[1:length(table(sceObject$sample))]
    picture5@left_annotation@anno_list[["sample"]]@fun@var_env[["color_mapping"]]@full_col[1:length(table(sceObject$sample))] <- samplecolor[1:length(table(sceObject$sample))]
    
    picture5@left_annotation@anno_list[["group"]]@color_mapping@full_col[1:length(table(sceObject$group))] <- color[1:length(table(sceObject$group))]
    picture5@left_annotation@anno_list[["group"]]@fun@var_env[["color_mapping"]]@full_col[1:length(table(sceObject$group))] <- color[1:length(table(sceObject$group))]
    draw(picture5)
    
    plotsave()
    
    plotfile(savepath = saveplotpath,
             mapname = "sce_scale_ExprHeatmap_no_clust.pdf",
             width = 10,height = 8,...)
    
    picture6 <- plotExprHeatmap(sceObject, scale = "last", row_anno = c("sample", "group"), row_clust = F,col_clust = F,
                                hm_pal = rev(hcl.colors(10, "YlGnBu")))
    picture6@left_annotation@anno_list[["sample"]]@color_mapping@full_col[1:length(table(sceObject$sample))] <- samplecolor[1:length(table(sceObject$sample))]
    picture6@left_annotation@anno_list[["sample"]]@fun@var_env[["color_mapping"]]@full_col[1:length(table(sceObject$sample))] <- samplecolor[1:length(table(sceObject$sample))]
    
    picture6@left_annotation@anno_list[["group"]]@color_mapping@full_col[1:length(table(sceObject$group))] <- color[1:length(table(sceObject$group))]
    picture6@left_annotation@anno_list[["group"]]@fun@var_env[["color_mapping"]]@full_col[1:length(table(sceObject$group))] <- color[1:length(table(sceObject$group))]
    draw(picture6)
    plotsave()
  
    df5 <- data.frame(picture5@matrix)
    rownames(df5) <- rownames(picture5@matrix)
    colnames(df5) <- colnames(picture5@matrix)
    
    savexlsx3(data = df5,
              filename = paste0(saveplotpath, "/heatmap_info.xlsx"),
              sheet = "heatmap_info")
  }

}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){

  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-sce","--sceObject",default = "rawdata/flowSet_object_sce.rds", help = "读取sce对象的rds文件名")
  parser$add_argument("-color","--color_by",default = "group", help = "绘图的颜色以该参数区分，可选：sample_id、condition、patient_id 默认group")
  parser$add_argument("-facet","--facet_col",default = 5, type = "integer", help = "密度图分面的列数")
  parser$add_argument("-group","--group_by",default = "sample", help = "确定绘图的分组，可选：sample_id、condition、patient_id 默认sample")
  parser$add_argument("-lab","--label_by",default = "sample", help = "定义图中的标签显示 默认sample")
  parser$add_argument("-scale","--scale_method",default = "last", help = "中位数标准化的方法，可选first、last、never 无指定,默认last")
  parser$add_argument("-fe","--features",default = NULL, help = "抗体类型，features = 'type' 表示谱系标签/ 'state' 表示状态标签/ NULL 无指定,默认NULL")
  
  parser$add_argument("-sp","--saveplotpath",default = "1.QC/", help = "保存图片的目录")
  parser$add_argument("-i","--dpi", default = 300, type = "integer", help = "图片清晰度")
  parser$add_argument("-fa","--family", default = "sans", help = "图片字体")
  
  args <- parser$parse_args()
  
  qcargs <- do.call(what = sce_QC_plot, args = args)
}
## 注意：每个目录和文件必须严格存在！rawdata是固定的下机数据存放目录 flowSet_object.rds是固定的文件名
