#!/opt/conda/bin/Rscript

#' 根据变量进行多元统计Permutation图可视化
#' 
#' @param prem 相应排序检测数据
#' @param mapname 图片名称
#' @param height 图片长度
#' @param width 图片宽度
#' @param other 其他参数
#' @param ... 见`ggplotsave`
#'
#' @export
permutationmap <- function(data,
                           mapname = "Permutation",
                           other = NULL,
                           width = 8, 
                           height = 5,
                           color = c("#446a37", "#19325f"),
                           fill = c("#c0d09d", "#106898"),
                           shape = c(21, 22),
                           ...) {

  suppressMessages(library("ggplot2"))
  options(warn=-1)
  
  if(is.null(data$prem)){
    return()
  }
  
  prem <- data$prem
  pp <- ggplot(data = prem$data, 
               mapping = aes(x = ID, y = va)) +
    geom_point(size = 3, 
               aes(shape = class, 
                   colour = class, 
                   fill = class)) +
    scale_shape_manual(values = shape) +
    scale_colour_manual(values = color) +
    scale_fill_manual(values = fill) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_abline(slope = prem$a1,
                intercept = prem$b1,
                linetype = "dashed") +
    geom_abline(slope = prem$a2,
                intercept = prem$b2,
                linetype = "dashed") +
    labs(title = prem$text,
         x = paste0(dim(prem$data)[1] / 2 - 1, " permutations 1 components"),
         y = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 15),
          panel.border = element_rect(colour = "black",linewidth = 1),
          # axis.line = element_line(colour = "black"),
          aspect.ratio = 3/4,
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"))
  
  args <- ggplotsave(plot = pp,
                     mapname = mapname,
                     other = other,
                     width = width,
                     height = height,
                     ...)
  return(args)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_mulstatistics_permutation <- map_autodraw$new(permutationmap)$draw

  parser <- ArgumentParser()
  parser$add_argument("-f","--filename", nargs="+",
                      help = "中间过程数据位置",required=T)
  parser$add_argument("-mn","--mapname",default = "Permutation", help = "保存名称")
  parser$add_argument("-sh","--shape",default = c(21,22),type= "integer", help = "点的形状",nargs=2)
  parser$add_argument("-col","--color",default = c("#446a37", "#19325f"), help = "点的边框颜色",nargs=2)
  parser$add_argument("-fl","--fill",default = c("#c0d09d", "#106898"), help = "点的填充颜色",nargs=2)
  
  # 基本参数
  parser$add_argument("-s","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf","html"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf","html"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  args <- parser$parse_args()

  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  flowargs <- do.call(what = map_mulstatistics_permutation,args = args)
}

#' 根据文件进行多元统计Permutation图可视化
#' 
#' @export
map_mulstatistics_permutation <- map_autodraw$new(permutationmap)$draw
