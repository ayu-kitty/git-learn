#!/opt/conda/bin/Rscript

#' 脂质组学特殊热图
#'
#' @param data 数据
#' @param legend 图例是否显示
#' @param min 热图数据最小值
#' @param max 热图数据最大值
#' @param breaks 热图断点范围
#' @param legend_breaks 图例是否断点标注
#' @param legend_labels 图例断点标注的标题
#' @param cellheight 热图单元格高度
#' @param cellwidth 热图单元格宽度
#' @param show_rownames 是否显示行名
#' @param show_colnames 是否显示列名
#' @param border_color 单元格边框颜色 默认无边框
#' @param fontsize_row 行名字体大小
#' @param fontsize_col 列名字体大小
#' @param angle_col 列名角度
#' @param display_numbers 是否显示文字以及显示什么文字
#' @param number_color 单元格文字颜色
#' @param number_format 单元格中数值的显示方式
#' @param fontsize_number 单元格内文字大小
#' @param na_col NA值对应的单元格填充颜色
#' @param color 热图颜色
#' @param mode 绘图类型 Carbon/Unsaturantion
#' @param mapname 图片保存名称
#' @param imagetype 图片保存格式
#' @param width 图片保存宽度
#' @param height 图片保存高度
#' @param family 字体,默认为Aria
#' @param dpi 分辨率
#' @param savepath 图片保存路径
#' 
#' @export
Lipid_pheatmap<-function(data,
                         legend = TRUE,
                         min=NULL,
                         max=NULL,
                         breaks=NULL,
                         legend_breaks = NULL,
                         legend_labels = NULL, 
                         cellheight = NULL,
                         cellwidth = NULL,
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         border_color = F,
                         fontsize_row =  NULL,
                         fontsize_col = NULL,
                         angle_col = 0,
                         display_numbers = NULL,
                         number_color= "black",
                         number_format = "%.2f",
                         fontsize_number =NULL,
                         na_col = "white",
                         color = colorRampPalette(c("navy","#366dbf", "#f1f7f3", "#e6412c","firebrick3"))(1000),
                         mode = "Carbon",
                         mapname=NULL,
                         imagetype=c("png","pdf"),
                         width = NULL,
                         height = NULL,
                         family = "sans",
                         dpi = 300,
                         savepath = "./",
                         ...){
  
  mat <- data
  
  if(is.null(cellheight)){ cellheight=ifelse(1200 / dim(mat)[1] > 15, 15, 1200 / dim(mat)[1] )}
  
  if(is.null(cellwidth)){  
    if(mode=="Carbon"){cellwidth=cellheight
    }else{cellwidth=2*cellheight}
  }
  
  if(is.null(fontsize_row)){ fontsize_row=0.75 * cellheight }
  
  if(is.null(fontsize_col)){ 
    if(mode=="Carbon"){fontsize_col=ifelse(0.5 * cellwidth > 15, 15, 0.75 * cellwidth) 
    }else{fontsize_col=ifelse(0.5 * cellwidth > 15, 15, 0.5 * cellwidth)}
  }
  
  if(is.null(width)){ 
    if(mode=="Carbon"){width=(dim(mat)[2] * cellwidth  + 200)/66
    }else{width=(dim(mat)[2] * cellwidth  + 200)/55}
  }
  
  if(is.null(height)){ height=(dim(mat)[1] * cellheight + 200)/72 }
  
  if(is.null(breaks)){ breaks= c(seq(min, 0, length.out=ceiling(500) + 1), seq(max/1000, max, length.out=floor(500)))}
  
  if(is.null(fontsize_number)){ fontsize_number=0.8*fontsize_col }
  
  plotfile(savepath = savepath,
           mapname = mapname,
           imagetype = imagetype,
           height = height,
           width = width,
           dpi = dpi,
           family = family,
           units = "in")
  
  pheatmap::pheatmap(mat = mat,
                     scale="none",
                     cluster_rows = F,
                     cluster_cols = F ,
                     legend = legend,
                     legend_breaks = legend_breaks,
                     legend_labels = legend_labels, 
                     cellheight = cellheight,
                     cellwidth = cellwidth,
                     show_rownames = show_rownames,
                     show_colnames = show_colnames,
                     main = NA,
                     border_color = border_color,
                     fontsize_row =  fontsize_row,
                     fontsize_col = fontsize_col,
                     angle_col = angle_col,
                     display_numbers = display_numbers,
                     number_color=number_color,
                     number_format = number_format,
                     fontsize_number = fontsize_number,
                     width = width,
                     height = height,
                     na_col = na_col,
                     color = color,
                     breaks=breaks,
                     ... )	
  
  plotsave()
  
  return(list(height = height,width = width))
} 			

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_lipid_heatmap <- map_autodraw$new(moudle = Lipid_pheatmap,row.names = 1)$draw
  
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
  
  #此图参数
  # parser$add_argument("-l","--legend",default = TRUE, help = "图例是否显示")
  # parser$add_argument("-min","--min",default = NULL, help = "热图数据最小值")
  # parser$add_argument("-max","--max",default = NULL, help = "热图数据最大值")
  # parser$add_argument("-br","--breaks",default = NULL, help = "热图断点范围")
  # parser$add_argument("-lb","--legend_breaks",default = NULL, help = "图例是否断点标注")
  # parser$add_argument("-ll","--legend_labels",default = NULL, help = "图例断点标注的标题")
  # parser$add_argument("-ch","--cellheight",default = NULL, help = "热图单元格高度")
  # parser$add_argument("-cw","--cellwidh",default = NULL, help = "热图单元格宽度")
  # parser$add_argument("-sr","--show_rownames",default = TRUE, help = "是否显示行名")
  # parser$add_argument("-sc","--show_colnames",default = TRUE, help = "是否显示列名")
  # parser$add_argument("-bc","--border_color",default = F, help = "单元格边框颜色 默认无边框")
  # parser$add_argument("-fr","--fontsize_row",default = NULL, help = "行名字体大小")
  # parser$add_argument("-fc","--fontsize_col",default = NULL, help = "列名字体大小")
  # parser$add_argument("-ag","--angle_col",default = 0, help = "列名角度")
  # parser$add_argument("-dn","--display_numbers",default = NULL, help = "是否显示文字以及显示什么文字")
  # parser$add_argument("-nc","--number_color",default = "black", help = "单元格文字颜色")
  # parser$add_argument("-nf","--number_format",default = "%.2f", help = "单元格中数值的显示方式")
  # parser$add_argument("-fnum","-- fontsize_number",default = NULL, help = "单元格内文字大小")
  # parser$add_argument("-nac","--na_col",default = "white", help = "NA值对应的单元格填充颜色")
  # parser$add_argument("-c","--color",default = colorRampPalette(c("navy","#366dbf", "#f1f7f3", "#e6412c","firebrick3"))(1000), help = "热图颜色")
  # parser$add_argument("-br","--breaks",default =NULL, help = "热图图例截断")
  # parser$add_argument("-mo","--mode",default ="Carbon", help = "绘图类型 Carbon/unsaturantion")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  result <- do.call(what = map_lipid_heatmap,args = args) 
}

#' @export
map_lipid_heatmap <- map_autodraw$new(moudle = Lipid_pheatmap,row.names = 1)$draw
