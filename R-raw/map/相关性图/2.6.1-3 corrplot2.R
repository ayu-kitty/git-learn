#!/opt/conda/bin/Rscript

#' 相关性图可视化
#'
#' @param height 图片长度
#' @param width 图片宽度
#' @param family 图片字体
#' @param savepath 保存路径
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param height 图片长度
#' @param width 图片宽度
#' @param dpi 图片分辨率
#' @param family 字体
#' @param ... 见corrplot::corrplot,以下均为corrplot::corrplot参数
#'
#' @export
corrplot2 <- function(corr,
                      p.mat = NULL,
                      savepath = "./",
                      mapname = "corrplot",
                      imagetype = c("jpg","pdf"),
                      dpi = 300,
                      width = 10,
                      height = 10,
                      family = "sans",
                      units = "in",
                      title = "\n Correlation",
                      tl.col = "black",
                      sig.level = c("", "", ""),
                      insig = "label_sig",
                      upmethod = "circle",
                      downmethod = "number",
                      col = SelectColors(palette = "navyfire",n = 100),
                      tl.cex = 10 / (4 + (dim(corr)[1]) / 8 + (max(nchar(row.names(corr)))) / 8),
                      order = ifelse(all(corr[row(corr) == col(corr)] == 1), "hclust", "original"),
                      pch.cex = 6 / (4 + (dim(corr)[1]) / 8 + (max(nchar(row.names(corr)))) / 8),
                      number.cex = 6 / (4 + (dim(corr)[1]) / 8 + (max(nchar(row.names(corr)))) / 8),
                      ...) {
  
  
  plotfile(savepath = savepath,
           mapname = mapname,
           imagetype = imagetype,
           height = height,
           width = width,
           dpi = dpi,
           family = family,
           units = units)
  
  data2 <- corrplot::corrplot(corr = corr,
                              type = "upper",
                              order = order,
                              title = title,
                              tl.col = tl.col,
                              tl.cex = tl.cex,
                              col = col,
                              diag = T,
                              p.mat = p.mat,
                              insig = insig,
                              pch.cex = pch.cex,
                              number.cex = number.cex,
                              sig.level = sig.level,
                              cl.pos = "r",
                              tl.pos = "lt",
                              method = upmethod,
                              ...)
  data2 <- corrplot::corrplot(corr = corr,
                              type = "lower",
                              order = order,
                              title = title,
                              tl.col = tl.col,
                              tl.cex = tl.cex,
                              col = col,
                              diag = F,
                              p.mat = p.mat,
                              insig = insig,
                              pch.cex = pch.cex,
                              number.cex = number.cex,
                              sig.level = sig.level,
                              add = T,
                              cl.pos = "n",
                              tl.pos = "n",
                              method = downmethod,
                              ...)
  
  
  plotsave()
  
  return(data2)
}

#' @export
map_common_corrplot2 <- function(mapmoudle = corrplot2,
                                 ...){
  
  map_autodraw$new(moudle = datatocorrplot,row.names = 1)$draw(mapmoudle = mapmoudle,
                                                               ...)
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "x表达矩阵文件",required = T)
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
  
  # 此图参数
  parser$add_argument("-c","--col",default = "navyfire",help = "配色设置")
  parser$add_argument("-m","--mode",default = "mode1", help = "显著性呈现方式选择,可选mode1,mode2,mode3,mode3,其他",
                      choices = c("mode1","mode2","mode3"))
  parser$add_argument("-um","--upmethod",default = "circle", help = "相关性点形状，可选circle,square,ellipse,number,shade,color,pie",
                      choices = c("circle","square","ellipse","number","shade","color","pie"))
  parser$add_argument("-dm","--downmethod",default = "number", help = "相关性点形状，可选circle,square,ellipse,number,shade,color,pie",
                      choices = c("circle","square","ellipse","number","shade","color","pie"))
  parser$add_argument("-ntx","--ntransx",default = T, help = "x数据是否转置",action='store_false',dest = "transx")
  parser$add_argument("-tx","--transx",default = T, help = "x数据是否转置",action='store_false')
  parser$add_argument("-pa","--adjust",default = "none", help = "显著性校正方法,默认none,包括holm,hochberg,hommel,bonferroni,BH,BY,fdr",
                      choices = c("none","holm","hochberg","hommel","bonferroni","BH","BY","fdr"))
  parser$add_argument("-cm","--cormethod",default = "pearson", help = "相关性计算方法,包括pearson,spearman,kendall",
                      choices = c("pearson","spearman","kendall"))
  parser$add_argument("-u","--use",default = "pairwise", help = "相关性计算方法,包括pairwise,complete",
                      choices = c("pairwise","complete"))
  
  args <- parser$parse_args()
  
  args$col <- SelectColors(palette = args$col,n = 100)
  
  args$insig <- switch(args$mode,
                       mode1 = "label_sig",
                       mode2 = "label_sig",
                       mode3 = "blank",
                       mode4 = "label_sig")
  args$sig.level <- switch(args$mode,
                           mode1 = c("", "", ""),
                           mode2 = c(0.01,0.05),
                           mode3 = 0.05,
                           mode4 = c(0.001,0.01,0.05))
  args$mode <- NULL
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  corrdata <- do.call(what = map_common_corrplot2,args = args)
  
}