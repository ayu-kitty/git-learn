#!/opt/conda/bin/Rscript

#' 空代质谱聚类图
#'
#' @param filename 文件路径
#' @param savepath 保存路径
#' @param mapname 保存名称
#' @param imagetype 图片格式
#' @param width 图片宽度
#' @param height 图片长度
#' @param asp 图像长宽比
#' @param lightmode 成像模式
#' @param xlab x轴标签
#' @param ylab y轴标签
#' @param family 字体
#' @param layout 分面分布
#' @param col 聚类颜色
#' @param lengedallshow 逻辑，聚类图例全部显示
#' @param mirrorx 逻辑，x轴镜像对称
#' @param mirrory 逻辑，y轴镜像对称
#' @param trans 逻辑，旋转90
#' @param ... 见[Cardinal::image()]
#'
#' @export
imzmlclusterimage <- function(filename,
                              savepath = "./",
                              mapname = "test",
                              imagetype = c("png","pdf"),
                              width = 8,
                              height = 8,
                              asp = 1,
                              lightmode = T,
                              xlab = "",
                              ylab = "",
                              family = "sans",
                              layout = NULL,
                              col = c("green4", "blue3", "firebrick",
                                      "gold", "darkviolet", "darkorange",
                                      "skyblue3", "olivedrab3", "dodgerblue3",
                                      "aquamarine2", "deeppink3", "slateblue3",
                                      "brown2", "palegreen2", "chocolate2",
                                      "antiquewhite3", "steelblue1", "violetred1",
                                      "burlywood3", "pink1", "slategray2",
                                      "orangered1", "cyan3", "yellow4",
                                      "red", "plum", "greenyellow",
                                      "mediumpurple2", "tan1", "magenta"),
                              lengedallshow = T,
                              mirrorx = NA,
                              mirrory = NA,
                              trans = NA,
                              values = "class",
                              figdata = NULL,
                              figfile = NULL,
                              ...) {
  print("质谱成像绘图开始")
  
  suppressMessages(library("Cardinal"))
  
  data <- readdata(filename = filename)

  ssc <- data$sscc
  
  samplename <- levels(run(ssc))
  
  if (is.null(layout)) {
    n1 <- ceiling(sqrt(length(samplename)))
    n2 <- ceiling(length(samplename) / n1)
    layout <- c(n2, n1)
  }else{
    n2 <- layout[1]
    n1 <- layout[2]
  }

  xmin <- min(coord(ssc)$x)
  ymin <- min(coord(ssc)$y)
  xlength1 <- NULL
  ylength1 <- NULL
  
  for ( i in 1:length(samplename)) {
    
    # i <- 1
    if(mirrorx[i] & !is.na(mirrorx[i])){
      xdata <- coord(ssc)$x[run(ssc) == samplename[i]]
      coord(ssc)$x[run(ssc) == samplename[i]] <- max(xdata)+min(xdata)-xdata
    }
    if(mirrory[i] & !is.na(mirrory[i])){
      ydata <- coord(ssc)$y[run(ssc) == samplename[i]]
      coord(ssc)$y[run(ssc) == samplename[i]] <- max(ydata)+min(ydata)-ydata
    }
    if(trans[i] & !is.na(trans[i])){
      coordxy <- coord(ssc)$x[run(ssc) == samplename[i]]
      coord(ssc)$x[run(ssc) == samplename[i]] <- coord(ssc)$y[run(ssc) == samplename[i]]
      coord(ssc)$y[run(ssc) == samplename[i]] <- coordxy
    }
    
    xlength1 <- c(xlength1,max(coord(ssc)$x[run(ssc) == samplename[i]])-min(coord(ssc)$x[run(ssc) == samplename[i]]))
    ylength1 <- c(ylength1,max(coord(ssc)$y[run(ssc) == samplename[i]])-min(coord(ssc)$y[run(ssc) == samplename[i]]))
    
  }
  
  xlength <- max(xlength1)
  xmax <- xmin+xlength
  ylength <- max(ylength1)
  ymax <- ymin+ylength
  
  plotfile(savepath = savepath,
           mapname = mapname,
           imagetype = imagetype,
           width = width * n1, height = height * n2,
           family = family,
           units = "px")
  
  if (lightmode) {
    lightmode()
    par(bg = "white")
  } else {
    darkmode()
  }
  showtext::showtext_auto()
  
  if(lengedallshow){
    
    print(image(ssc,
                values = values,
                xlab = xlab,
                ylab = ylab,
                xlim = c(xmin - xlength * 0.05, xmax + xlength * 0.2),
                ylim = c(ymin - ylength * 0.1, ymax + ylength * 0.1),
                asp = asp,
                layout = layout,
                col = col,
                ...))
    
    
  }else{
    sampleclass <- ssc$class[[1]]
    sampleclassname <- levels(run(ssc))
    i <- 1
    samplerange <<- (run(ssc) == sampleclassname[i])
    samplerangeclass <- unique(sampleclass[run(ssc) == sampleclassname[i]])
    samplerangeclass <- samplerangeclass[order(samplerangeclass)]
    col1 <- col[samplerangeclass]
    key <- list(legend = samplerangeclass,fill = col1)
    print(image(ssc,
                values = "class",
                xlab = xlab,
                ylab = ylab,
                xlim = c(xmin - xlength * 0.05, xmax + xlength * 0.2),
                ylim = c(ymin - ylength * 0.1, ymax + ylength * 0.1),
                asp = asp,
                layout = layout,
                col = col,
                subset = samplerange,
                key = key,
                ...))
    if(length(sampleclassname) > 1){
      for ( i in 2:length(sampleclassname)) {
        samplerange <<- (run(ssc) == sampleclassname[i])
        samplerangeclass <- unique(sampleclass[run(ssc) == sampleclassname[i]])
        samplerangeclass <- samplerangeclass[order(samplerangeclass)]
        col1 <- col[samplerangeclass]
        key <- list(legend = samplerangeclass,fill = col1)
        print(image(ssc,
                    values = "class",
                    xlab = xlab,
                    ylab = ylab,
                    xlim = c(xmin - xlength * 0.05, xmax + xlength * 0.2),
                    ylim = c(ymin - ylength * 0.1, ymax + ylength * 0.1),
                    asp = asp,
                    layout = NULL,
                    col = col,
                    subset = samplerange,
                    key = key,
                    ...))
      }
    }
  }
  
  if(!is.null(figfile)){
    figdata <- png::readPNG(figfile)
  }
  if(!is.null(figdata)){
    rasterImage(figdata,
                xleft = xmin,
                ybottom = ymax,
                xright = xmax,
                ytop = ymin)
  }
  
  plotsave()
  showtext::showtext_auto(FALSE)
  
  gc(reset = TRUE)
  return("完成绘图")
}


#' 空代聚类statistic图
#'
#' @param filename 文件路径
#' @param savepath 保存路径
#' @param mapname 保存名称
#' @param imagetype 图片格式
#' @param width 图片宽度
#' @param height 图片长度
#' @param lightmode 成像模式
#' @param family 字体
#' @param col 聚类颜色
#' @param ... 见[Cardinal::plot()]
#' @param values 
#' @param israwdata 
#'
#' @export
imzmlclusterplot <- function(filename,
                             savepath = "./",
                             mapname = "test",
                             imagetype = c("jpg","pdf"),
                             width = 16,
                             height = 9,
                             lightmode = T,
                             family = "sans",
                             col = c("green4", "blue3", "firebrick",
                                     "gold", "darkviolet", "darkorange",
                                     "skyblue3", "olivedrab3", "dodgerblue3",
                                     "aquamarine2", "deeppink3", "slateblue3",
                                     "brown2", "palegreen2", "chocolate2",
                                     "antiquewhite3", "steelblue1", "violetred1",
                                     "burlywood3", "pink1", "slategray2",
                                     "orangered1", "cyan3", "yellow4",
                                     "red", "plum", "greenyellow",
                                     "mediumpurple2", "tan1", "magenta"),
                             values = "statistic",
                             ...) {
  print("质谱成像绘图开始")
  
  suppressMessages(library("Cardinal"))
  
  data <- readRDS(file = filename)
  
  ssc <- data$sscc
  
  plotfile(savepath = savepath,
           mapname = mapname,
           imagetype = imagetype,
           width = width, height = height,
           family = family,
           units = "px")
  
  if (lightmode) {
    lightmode()
    par(bg = "white")
  } else {
    darkmode()
  }
  showtext::showtext_auto()
  print(
    Cardinal::plot(ssc,
                   values = values,
                   col = col,
                   ...
    )
  )
  
  plotsave()
  showtext::showtext_auto(FALSE)
  
  gc(reset = TRUE)
  return("完成绘图")
}
