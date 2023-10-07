#!/opt/conda/bin/Rscript

#' 根据变量进行多元统计3Dscore图可视化
#' 
#' @param data 数据
#' @param mapname 图片名称
#' @param height 图片长度
#' @param width 图片宽度
#' @param x x轴标签
#' @param y y轴标签
#' @param z z轴标签
#' @param title 标题
#' @param fill 点填充颜色
#' @param shape 点形状
#' @param family 字体
#' @param mode 多元统计模式
#' @param size 点大小
#' @param xlim x轴范围
#' @param ylim y轴范围
#' @param zlim z轴范围
#'
#' @export
mulstatistics3dmap <- function(data,
                               mapname = "statistics",
                               x = colnames(data)[2],
                               y = colnames(data)[3],
                               z = colnames(data)[4],
                               family = "sans",                              
                               mode = "PCA",
                               title = mode,
                               point_size = if(dim(data)[1] > 500){I(10)}else{I(200)},
                               xlim = c(-max(abs(data[, 2])) * 1.1, max(abs(data[, 2])) * 1.1),
                               ylim = c(-max(abs(data[, 3])) * 1.1, max(abs(data[, 3])) * 1.1),
                               zlim = c(-max(abs(data[, 4])) * 1.1, max(abs(data[, 4])) * 1.1),
                               classfile = "classtype.xlsx",
                               point_fill = stylefun_group(classfile = classfile,object = data,styletype = "fill"),
                               point_shape = stylefun_group(classfile = classfile,object = data,styletype = "shape"),
                               width = 8,
                               height = 6,
                               savepath = "./",
                               imagetype = "html",
                               ...) {
  
  suppressMessages(library("plotly"))
  options(warn=-1)
  
  score <- data
  names(score)[2] <- c("p1")
  names(score)[3] <- c("p2")
  names(score)[4] <- c("p3")
  
  if (is.numeric(point_shape)) {
    point_shape[point_shape == 24] <- 1
    point_shape[point_shape == 25] <- 0
  }
  
  score$Group <- factor(score$Group, levels = unique(score$Group))
  
  p <- plot_ly(score,
               type = "scatter3d",
               x = ~p1, y = ~p2, z = ~p3,
               color = ~Group, colors = point_fill,
               size = point_size,
               # stroke = ~group, strokes=colour,
               symbol = ~Group, symbols = point_shape,
               mode = "markers",
               hoverinfo = "text",
               text = paste(
                 "Sample:", score[,1],
                 "<br>Group:", score$Group,
                 "<br>PC1:", score$p1,
                 ifelse(mode == "OPLS", "<br>PCo1:", "<br>PC2:"), score$p2,
                 ifelse(mode == "OPLS", "<br>PCo2:", "<br>PC3:"), score$p3
               ),
               width = width*100, height = height*100) %>%
    plotly::layout(title = paste0("<br>",title),
                   font = list(family = family),
                   scene = list(xaxis = list(title = x,
                                             gridcolor = "#999999",
                                             range = xlim,
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwidth = 2),
                                yaxis = list(title = y,
                                             gridcolor = "#999999",
                                             range = ylim,
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2),
                                zaxis = list(title = z,
                                             gridcolor = "#999999",
                                             range = zlim,
                                             zerolinewidth = 1,
                                             ticklen = 5,
                                             gridwith = 2),
                                aspectratio = list(x = 0.9, y = 0.9, z = 0.9)
                   ),
                   paper_bgcolor = "#FFFFFF",
                   plot_bgcolor = "#FFFFFF")
  
  if(!is.null(imagetype)){
    for (type in imagetype) {
      if(is.na(type)){
        print(p)
      }else if(type == "html"){
        # runinpath(path = savepath,
        #           moudle = htmltools::save_html,
        #           moudlename = "html保存", 
        #           html = p,
        #           file = paste0(mapname,".html"))
        
        logfile <- file(tempfile(), open = "wt")
        sink(file = logfile)
        sink(file = logfile,type = "message")
        runinpath(path = savepath,
                  moudle = htmlwidgets::saveWidget,
                  moudlename = "html保存", 
                  widget = p,
                  file = paste0(mapname,".html"),
                  selfcontained = TRUE)
        sink(type = "message")
        sink()
      }else{
        loadpython()
        runinpath(path = savepath,
                  moudle = save_image,
                  moudlename = "图片保存", 
                  p = p,
                  file = paste0(mapname,".",type),
                  scale = 5)
      }
    }
  }
  
  return(list(plot=p))
}

#' @export
mulstatistics3dmap_obj <- function(data,...){
  if("mulstatistics" %in% class(data)){
    data <- getscore_file(rdspath = data)
  }
  result <- mulstatistics3dmap(data = data,...)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_mulstatistics_3dscoremap <- map_autodraw$new(mulstatistics3dmap_obj)$draw
  
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
  parser$add_argument("-ps","--point_shape",default = "", help = "点形状",nargs="+")
  parser$add_argument("-pi","--point_size",default = 200, type= "double", help = "点大小")
  
  parser$add_argument("-mo","--mode",default = "")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}

  if("" %in% args$point_fill){args$point_fill <- NULL}
  if("" %in% args$point_shape){args$point_shape <- NULL}
  args$point_size <- I(args$point_size)
  
  result <- do.call(what = map_mulstatistics_3dscoremap,args = args) 
}

#' 根据文件进行多元统计3Dscore图可视化
#' 
#' @export
map_mulstatistics_3dscoremap <- map_autodraw$new(mulstatistics3dmap_obj)$draw
