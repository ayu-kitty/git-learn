#!/opt/conda/bin/Rscript

#' corrnetwork
#'
#' 相关性网络图
#'
#' @param height 图片长度
#' @param width 图片宽度
#' @param family 图片字体
#' @param color 颜色
#' @param layout 点分布模式，见[ggraph()]中layout
#' @param algorithm 点分布模式，见[ggraph()]中algorithm
#' @param linewidth 线粗细范围
#' @param linealpha 线透明度
#' @param textsize 文本标签大小
#' @param nodealpha 点透明度
#' @param repel 逻辑，是否启用repel
#' @param fontface 标签吸引力
#' @param max.overlaps 名称重叠度
#' @param nodedata 点数据
#' @param nodesize 点大小范围
#' @param shape 点形状
#' @param linkdata 连接数据
#' @param imagetype 图片格式
#' @param plot.margin 边缘空白
#' @param ... 见ggplotsave
#'
#' @export
corrnetwork <- function(linkdata,
                        nodedata,
                        mapname = "corrnetwork",
                        imagetype = c("jpg","pdf"),
                        width = 15,
                        height = 15,
                        layout = "igraph",
                        algorithm = "circle",
                        color = SelectColors(palette = "navyfire",
                                             object = c("negative correlation","white","positive correlation"),
                                             subset = c(1,3)),
                        linewidth = c(0.2, 1),
                        nodesize = c(2,6),
                        linealpha = 0.3,
                        textsize = 2 + 30 / length(unique(c(linkdata[, 1], linkdata[, 2]))),
                        nodealpha = 0.6,
                        fontface = "bold",
                        shape = c(21,22),
                        fill = SelectColors(palette = "blindless",n = 2),
                        repel = T,
                        max.overlaps = 100,
                        plot.margin = margin(100,100,100,100),
                        family = "sans",
                        ...) {
  
  options(warn = -1)
  suppressMessages(library("ggraph"))
  suppressMessages(library("igraph"))

  g <- graph_from_data_frame(d = linkdata, 
                             directed = F,
                             vertices = nodedata)
  pp <- ggraph(g, layout = layout, algorithm = algorithm) +
    geom_edge_link(mapping = aes(edge_width = abs(cor), 
                                 edge_color = class),
                   alpha = linealpha,
                   show.legend = F) +
    geom_node_point(mapping = aes(shape = type,
                                  fill = type,
                                  size = Freq),
                    show.legend = F, 
                    alpha = nodealpha) +
    geom_node_text(mapping = aes(label = name),
                   repel = repel, size = textsize,
                   family = family,
                   fontface = fontface,
                   max.overlaps = max.overlaps) +
    theme_graph(plot_margin = plot.margin,
                base_family = family) +
    scale_edge_width_continuous(range = linewidth) +
    scale_edge_color_manual(values = color)+
    scale_shape_manual(values = shape)+
    scale_fill_manual(values = fill)+
    scale_size_continuous(range = nodesize)
  
  plotdata <- ggplotsave(plot = pp,
                         imagetype = imagetype,
                         width = width,
                         height = height,
                         family = family,
                         mapname = mapname,
                         ...)
  return(plotdata)
}

#' @export
datatocorrnetwork <- function(filename = "1",
                              data = readdata(filename = filename,row.names = 1),
                              filenamey = NULL,
                              ci = FALSE,
                              xname = "Featurex",
                              yname = "Featurey",
                              corfilter = 0.95,
                              corfiltertype = "+-",
                              pfilter = 0.05,
                              adjust = "none",
                              transx = T,
                              transy = T,
                              cormethod = "pearson",
                              use = "pairwise",
                              mapmoudle = corrnetwork,
                              mapname = "corrnetwork",
                              savepath = "./",
                              savecorr = T,
                              ...){
  
  x <- data
  y <- readdata(filename = filenamey,row.names = 1)
  
  if(is.null(data)){
    
    savetxt(data = "数据为空,不进行相关性分析",
            filename = paste0(savepath,"/说明.txt"))
    
    return()
    
  }
  
  if(dim(data)[1] < 3 & is.null(y)){
    
    savetxt(data = "数据太少,不进行相关性分析",
            filename = paste0(savepath,"/说明.txt"))
    
    return()
    
  }
  
  if(dim(data)[2] < 3){
    
    savetxt(data = "数据太少,不进行相关性分析",
            filename = paste0(savepath,"/说明.txt"))
    
    return()
    
  }
  
  corrdata <- corrcal(x = x,
                      y = y,
                      ci = ci,
                      xname = xname,
                      yname = yname,
                      corfilter = corfilter,
                      corfiltertype = corfiltertype,
                      pfilter = pfilter,
                      adjust = adjust,
                      transx = transx,
                      transy = transy,
                      method = cormethod,
                      use = use)
  
  if(is.null(corrdata)){
    savetxt(data = "数据为空,不进行相关性分析",
            filename = paste0(savepath,"/说明.txt"))
    
    return(corrdata)
  }else if(dim(corrdata$filterdata$linkdata)[1] == 0){
    
    savetxt(data = paste0("经过相关性>",corfilter,",显著性<",pfilter,",无符合条件的结果"),
            filename = paste0(savepath,"/说明.txt"))
    
    return(corrdata)
  }
  
  plotdata <- mapmoudle(linkdata = corrdata$filterdata$linkdata,
                        nodedata = corrdata$filterdata$nodedata,
                        mapname = mapname,
                        savepath = savepath,
                        ...)
  
  corrdata$plotdata <- plotdata
  
  if (savecorr) {
    savexlsx1(data = corrdata$filterdata$linkdata, 
              filename = paste0(savepath,"/",mapname,".xlsx"), 
              sheet = "cordata")
  } 
  
  return(corrdata)
  
}


#' @export
map_common_corrnetwork <- function(mapmoudle = corrnetwork,
                                   ...){
  
  map_autodraw$new(moudle = datatocorrnetwork,row.names = 1)$draw(mapmoudle = mapmoudle,
                                                                  ...)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "x表达矩阵文件",required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  parser$add_argument("-fy","--filenamey",default = NULL,nargs="*",
                      help = "y表达矩阵文件,如果是xlsx文件以`数据矩阵.xlsx 数据矩阵`形式传参")
  
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
  parser$add_argument("-co","--color",default = "navyfire",help = "配色设置")
  parser$add_argument("-xn","--xname",default = "Featurex", help = "x表格类型")
  parser$add_argument("-yn","--yname",default = "Featurey", help = "y表格类型")
  parser$add_argument("-cf","--corfilter",default = 0.95, type= "double",help = "相关性筛选标准")
  parser$add_argument("-pf","--pfilter",default = 0.05, type= "double",help = "显著性筛选标准")
  parser$add_argument("-l","--layout",default = "igraph",
                      help = "点分布模块,包括auto,igraph,dendrogram,manual,linear,matrix,treemap,circlepack,partition,hive")
  parser$add_argument("-a","--algorithm",default = "circle",
                      help = "点分布形状,包括tree,sugiyama,bipartite,star,circle,nicely,dh,gem,graphopt,grid,mds,sphere,randomly,fr,kk,drl,lgl等")
  parser$add_argument("-lw","--linewidth",default = c(0.2,1),nargs=2,type= "double",help = "线粗细范围")
  parser$add_argument("-ns","--nodesize",default = c(2,6),nargs=2,type= "double",help = "点大小范围")
  parser$add_argument("-ts","--textsize",default = 0,type= "double",help = "文本大小")
  parser$add_argument("-ntx","--ntransx",default = T, help = "x数据是否转置",action='store_false',dest = "transx")
  parser$add_argument("-nty","--ntransy",default = T, help = "y数据是否转置",action='store_false',dest = "transy")
  parser$add_argument("-pa","--adjust",default = "none", help = "显著性校正方法,默认none,包括holm,hochberg,hommel,bonferroni,BH,BY,fdr",
                      choices = c("none","holm","hochberg","hommel","bonferroni","BH","BY","fdr"))
  parser$add_argument("-cm","--cormethod",default = "pearson", help = "相关性计算方法,包括pearson,spearman,kendall",
                      choices = c("pearson","spearman","kendall"))
  parser$add_argument("-u","--use",default = "pairwise", help = "相关性计算方法,包括pairwise,complete",
                      choices = c("pairwise","complete"))
  
  args <- parser$parse_args()
  
  args$color <- SelectColors(palette = args$color,
                             object = c("positive correlation","white","negative correlation"),
                             subset = c(1,3))
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  if(args$textsize == 0){ args$textsize <- NULL}
  
  corrdata <- do.call(what = map_common_corrnetwork,args = args)
  
}
