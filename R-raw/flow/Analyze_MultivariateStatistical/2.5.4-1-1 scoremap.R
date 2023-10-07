#!/opt/conda/bin/Rscript

#' 根据变量进行多元统计score图可视化
#' 
#' @param data 数据
#' @param mapname 图片名称
#' @param height 图片长度
#' @param width 图片宽度
#' @param x x轴标签
#' @param y y轴标签
#' @param title 标题
#' @param showname 逻辑，是否显示样本名
#' @param max.overlaps 名称重叠度
#' @param fill 点填充颜色
#' @param colour 点边框颜色
#' @param shape 点形状
#' @param other 其他参数
#' @param ... 见`ggplotsave`
#'
#' @export
mulstatisticsmap <- function(data,
                             mapname = "statistics",
                             x = colnames(data)[2],
                             y = colnames(data)[3],
                             mode = "",
                             title = mode,
                             classfile = "classtype.xlsx",
                             point_fill = stylefun_group(classfile = classfile,object = data,styletype = "fill"),
                             point_colour = stylefun_group(classfile = classfile,object = data,styletype = "colour"),
                             point_shape = stylefun_group(classfile = classfile,object = data,styletype = "shape"),
                             point_size = 4,
                             showname = F,
                             max.overlaps = 200,
                             other = NULL,
                             width = 8,
                             height = 5,
                             imagetype = c("jpg","pdf","html"),
                             stat_ellipse = ggplot2::stat_ellipse(geom = "path",
                                                                  type = "norm",
                                                                  alpha = 1,
                                                                  segments = 200, 
                                                                  show.legend = F,
                                                                  szie = 0.1),
                             hline = ggplot2::geom_hline(yintercept = 0,linetype = 2,col = "grey"),
                             vline = ggplot2::geom_vline(xintercept = 0,linetype = 2,col = "grey"),
                             addfacet = F,
                             facetlistname = "Sample",
                             facetlayout = NULL,
                             shrink = T,
                             scales = "fixed",
                             gradientmode = F,
                             gradientshape = 21,
                             gradientalpha = 0.5,
                             gradientcolours = colorRampPalette(c("grey", "blue2", "cyan2", "yellow","brown1", "firebrick3"))(1000),
                             ...) {
  suppressMessages(library("ggplot2"))
  suppressMessages(library("ggforce"))
  options(warn=-1)
  
  score <- data
  names(score)[2] <- c("p1")
  names(score)[3] <- c("p2")
  
  if (gradientmode) {
  } else {
    score$Group <- factor(score$Group, levels = unique(score$Group))
  }
  
  if (addfacet) {
    score[, "facet"] <- score[, facetlistname]
  }
  
  if (gradientmode) {
    pp <- ggplot(data = score, mapping = aes(x = p1, 
                                             y = p2))+
      labs(title = title, x = x, y = y)+
      hline+
      vline+
      stat_ellipse+
      geom_point(size = point_size,
                 stroke = point_stroke,
                 shape = gradientshape,
                 alpha = gradientalpha,
                 mapping = aes(fill = Group,
                               colour = Group,
                               text = paste0("Sample:",Sample, "<br>Intensity:", Group)))+
      scale_fill_gradientn(colours = gradientcolours) +
      scale_colour_gradientn(colours = gradientcolours)+
      theme_bw()+
      theme(legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 15),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black",linewidth = 1),
            # axis.line = element_line(colour = "black"),
            aspect.ratio = 3/4)
  }else{
    pp <- ggplot(data = score, mapping = aes(x = p1, 
                                             y = p2))+
      labs(title = title, x = x, y = y)+
      hline+
      vline+
      stat_ellipse+
      geom_point(size = point_size,
                 mapping = aes(shape = Group,
                               fill = Group,
                               colour = Group,
                               text = paste0("Sample:",Sample, "<br>Group:", Group)))+
      scale_shape_manual(values = point_shape)+
      scale_fill_manual(values = point_fill)+
      scale_colour_manual(values = point_colour)+
      theme_bw()+
      theme(legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 15),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black",linewidth = 1),
            # axis.line = element_line(colour = "black"),
            aspect.ratio = 3/4)
  } 
  
  if (showname) {
    pp <- pp + ggrepel::geom_text_repel(mapping = aes(label = Sample),
                                        max.overlaps = max.overlaps,
                                        na.rm = T)
  }
  
  if (addfacet) {
    facetnumber <- length(unique(score[, "facet"]))
    if (facetnumber > 1) {
      if (is.null(facetlayout)) {
        n1 <- ceiling(sqrt(facetnumber))
        n2 <- ceiling(facetnumber / n1)
      } else {
        n1 <- facetlayout[2]
        n2 <- facetlayout[1]
      }
      pp <- pp + facet_wrap(~facet, scales = scales, shrink = shrink, nrow = n2, ncol = n1)
    } else {
      n1 <- 1
      n2 <- 1
    }
  } else {
    n1 <- 1
    n2 <- 1
  }
  
  args <- ggplotsave(plot = pp,
                     mapname = mapname,
                     other = other,
                     width = width*n1,
                     height = height*n2,
                     imagetype = imagetype,
                     ...)
  return(args)
}

#' @export
mulstatisticsmap_obj <- function(data,...){
  if("mulstatistics" %in% class(data)){
    data <- getscore_file(rdspath = data)
  }
  result <- mulstatisticsmap(data = data,...)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_mulstatistics_scoremap <- map_autodraw$new(mulstatisticsmap_obj)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "score矩阵文件",required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf","html"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf","html"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-c","--classfile",default = "classtype.xlsx", help = "模板")
  parser$add_argument("-fi","--point_fill",default = "", help = "填充颜色",nargs="+")
  parser$add_argument("-co","--point_colour",default = "", help = "点外圈颜色",nargs="+")
  parser$add_argument("-ps","--point_shape",default = "", help = "点形状",nargs="+")
  parser$add_argument("-pi","--point_size",default = 4, type= "double",help = "点大小")
  parser$add_argument("-sn","--showname",default = F, type="logical",help = "是否显示样本名")
  parser$add_argument("-mo","--max_overlaps",default = 200, type= "double",help = "样本名重叠度",dest="max.overlaps")
  parser$add_argument("-md","--mode",default = "")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  if("" %in% args$point_fill){args$point_fill <- NULL}
  if("" %in% args$point_colour){args$point_colour <- NULL}
  if("" %in% args$point_shape){args$point_shape <- NULL}
  
  result <- do.call(what = map_mulstatistics_scoremap,args = args)
}

#' 根据文件进行多元统计score图可视化
#' 
#' @export
map_mulstatistics_scoremap <- map_autodraw$new(mulstatisticsmap_obj)$draw
