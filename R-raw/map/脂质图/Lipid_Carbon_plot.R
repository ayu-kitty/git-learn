#!/opt/conda/bin/Rscript

#' @export
Lipid_Carbon_plot <- function(data,
                              savepath = "./",
                              mapname = "Lipid_Carbon",
                              wbsave = T,
                              savename = paste0(mapname,".xlsx"),
                              ...){
  
  suppressMessages(library("gtools"))
  suppressMessages(library("reshape2"))
  
  if("diff_filter" %in% class(data)){
    mat <- data
    plotdata <- merge(mat$args$info[,c("Sub Class","Total Carbon"),drop=F],
                      mat$args$data,by = 0,all = F)[,-1]
    data_C <- aggregate(.~`Sub Class`+`Total Carbon`,plotdata,sum) 
    
    diffdata <- fc_ana(data = data_C[,-1:-2,drop = F],
                       class = mat$args$singleclass,
                       compare = mat$args$compare,
                       log = mat$args$log)
    
    diffdata2 <- p_ana(data = data_C[,-1:-2,drop = F],
                       class = mat$args$singleclass,
                       paired = mat$args$paired,
                       p_method = mat$args$p_method,
                       p_adjust_methods = mat$args$p_adjust_methods)
    diffdata <- merge(diffdata,diffdata2,by = 0,sort = F,all = T)
    row.names(diffdata) <- diffdata[,1]
    diffdata <- diffdata[,-1]
    data_C <- merge(data_C,diffdata,by = 0,sort = F,all = T)[,-1]
    data_C[,"level"] <- ""
    data_C[data_C$`p-value` < 0.05,"level"] <- "*"
    data_C[data_C$`p-value`  < 0.01,"level"] <- "**"
    data_C[data_C$`p-value`  < 0.001,"level"] <- "***"
    data_C[data_C$`p-value`  < 0.0001,"level"] <- "****"
    result <- data_C
  }else{
    result <- data
  }
  
  resultp <- result[,c("Sub Class","Total Carbon","log2FoldChange")]
  resultlevel <- result[,c("Sub Class","Total Carbon","level")]
  plot <- dcast(resultp,`Sub Class`~`Total Carbon`,value.var = "log2FoldChange")
  plotlevel <- dcast(resultlevel,`Sub Class`~`Total Carbon`,value.var = "level")
  rownames(plotlevel) <- plotlevel[,1]
  plotlevel[is.na(plotlevel)] <- ""
  plotlevel <- plotlevel[,-1]
  rownames(plot) <- plot[,1]
  plot <- plot[,-1]
  plot_na <- plot
  plot_na[is.na(plot_na)] <- 0
  min <- min(plot_na)
  max <- max(plot_na)
  
  Lipid_pheatmap(data = plot,
                 mapname = mapname,
                 min = min,max=max,
                 display_numbers = plotlevel,
                 mode = "Carbon",
                 savepath = savepath,
                 ...)
  
  if(wbsave){
    savexlsx1(data = result,
              filename = paste0(savepath,"/",savename))
  }			  
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_lipid_carbonheatmap <- map_autodraw$new(moudle = Lipid_Carbon_plot)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "矩阵文件",
                      required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # parser$add_argument("-l","--legend",default = TRUE, help = "图例是否显示")
  # parser$add_argument("-lb","--legend_breaks",default = NULL, help = "图例是否断点标注")
  # parser$add_argument("-ll","--legend_labels",default = NULL, help = "图例断点标注的标题")
  # parser$add_argument("-sr","--show_rownames",default = TRUE, help = "是否显示行名")
  # parser$add_argument("-sc","--show_colnames",default = TRUE, help = "是否显示列名")
  # parser$add_argument("-ch","--cellheight",default = NULL, help = "热图单元格高度")
  # parser$add_argument("-cw","--cellwidh",default = NULL, help = "热图单元格宽度")
  # parser$add_argument("-bc","--border_color",default = F, help = "单元格边框颜色 默认无边框")
  # parser$add_argument("-fnum","-- fontsize_number",default = NULL, help = "单元格内文字大小")
  # parser$add_argument("-fr","--fontsize_row",default = NULL, help = "行名字体大小")
  # parser$add_argument("-fc","--fontsize_col",default = NULL, help = "列名字体大小")
  # parser$add_argument("-w","--width",default = NULL, help = "保存图片宽度")
  # parser$add_argument("-he","--height",default = NULL, help = "保存图片高度")
  # parser$add_argument("-c","--color",default = NULL, help = "热图颜色")
  # parser$add_argument("-dpi","--dpi",default = 300, help = "分辨率")
  # parser$add_argument("-ag","--angle_col",default = 0, help = "列名角度")
  # parser$add_argument("-nc","--number_color",default = "black", help = "单元格文字颜色")
  # parser$add_argument("-nf","--number_format",default = "%.2f", help = "单元格中数值的显示方式")
  # parser$add_argument("-tp","--type",default = "jpg,pdf", help = "图片保存格式")
  # parser$add_argument("-nac","--na_col",default = "white", help = "NA值对应的单元格填充颜色")
  # parser$add_argument("-p","--paletteLength",default = 1000, help = "取色长度")
  # parser$add_argument("-fa","--family",default = "sans", help = "字体")
  # parser$add_argument("-sp","--savepath",default = "./", help = "图片保存路径")
  # parser$add_argument("-fp","--frompath",default = "./", help = "数据来源路径")
  # parser$add_argument("-sn","--savename",default = "绘图数据-Carbon.xlsx", help = "绘图数据保存名称")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  result <- do.call(what = map_lipid_carbonheatmap,args = args) 
}  

#' @export
map_lipid_carbonheatmap <- map_autodraw$new(moudle = Lipid_Carbon_plot)$draw
