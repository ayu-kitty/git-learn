#!/opt/conda/bin/Rscript

#' 空代质谱成像图
#'
#' @param filename 文件路径
#' @param savepath 保存路径
#' @param mapname 保存名称
#' @param imagetype 图片格式
#' @param width 图片宽度
#' @param height 图片长度
#' @param asp 图像长宽比
#' @param lightmode 成像模式
#' @param fun 处理函数
#' @param smooth.image 是否平滑处理
#' @param superpose 是否分面绘制
#' @param normalize.image 归一化模式
#' @param contrast.enhance 成像强对比
#' @param xlab x轴标签
#' @param ylab y轴标签
#' @param color 成像颜色
#' @param family 字体
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param subset 区域选择
#' @param vague 是否虚化
#' @param iteration 虚化次数
#' @param maxPixels 虚化上限
#' @param area 是否绘制选择区域
#' @param areards 选择区域信息路径
#' @param areaname 选择区域名称
#' @param mirrorx 逻辑，x轴镜像对称
#' @param mirrory 逻辑，y轴镜像对称
#' @param trans 逻辑，旋转90
#' @param layout 多样本排布方式
#' @param intsenityrange 标签范围，c(0,10000)
#' @param intsenityratiorange 标签比例范围，c(0,0.2)
#' @param mapmz 逻辑值,是否绘制mz的图
#' @param attach.only 逻辑，数据是否读取到内存
#' @param zlim 图例范围
#' @param figdata HE图片数据
#' @param figfile HE图片路径
#' @param addy 是否按照y轴添加
#' @param addx 是否按照x轴添加
#' @param filtermz 进行mz筛选
#' @param ... 见[Cardinal::image()]
#'
#' @export
imzmlimage <- function(filename,
                       savepath = "./",
                       mapname = "test",
                       imagetype = c("png","pdf"),
                       width = 8,
                       height = 7,
                       asp = 1,
                       lightmode = F,
                       fun = mean,
                       smooth.image = "gaussian",
                       superpose = F,
                       normalize.image = "none",
                       contrast.enhance = "none",
                       xlab = "",
                       ylab = "",
                       color = c("blue2", "cyan2", "yellow",
                                 "brown1", "firebrick3"),
                       family = "sans",
                       mass.range = NULL,
                       resolution = 5,
                       units = "ppm",
                       subset = T,
                       vague = F,
                       iteration = 2,
                       maxPixels = 1000000,
                       area = F,
                       areards = NA,
                       areaname = gsub(pattern = "\\.rds$",
                                       replacement = "",
                                       basename(areards)),
                       intsenityrange = NULL,
                       intsenityratiorange = NULL,
                       zlim = NULL,
                       mapmz = F,
                       mz = NULL,
                       col = "black",
                       mirrorx = rep(F, ifelse(is.vector(filename),length(filename),1)),
                       mirrory = rep(F, ifelse(is.vector(filename),length(filename),1)),
                       trans = rep(F, ifelse(is.vector(filename),length(filename),1)),
                       layout = NULL,
                       attach.only = ifelse(is.null(intsenityratiorange) & is.null(intsenityrange), T, F),
                       figdata = NULL,
                       figfile = NULL,
                       addy = F,
                       addx = F,
                       filtermz = NULL,
                       ...) {
  print("成像图绘图开始")

  suppressMessages(library("Cardinal"))
  setCardinalBPPARAM(MulticoreParam(workers = 3))
  
  mse <- readdata(filename = filename,
                  mass.range = mass.range,
                  resolution = resolution,
                  units = units,
                  vague = vague,
                  iteration = iteration,
                  maxPixels = maxPixels,
                  area = area,
                  areards = areards,
                  areaname = areaname,
                  mirrorx = mirrorx,
                  mirrory = mirrory,
                  trans = trans,
                  addy = addy,
                  addx = addx,
                  filtermz = filtermz,
                  intsenityrange = intsenityrange,
                  intsenityratiorange = intsenityratiorange,
                  smooth.image = smooth.image,
                  attach.only = attach.only)
  
  if(vague){
    smooth.image <- "none"
  }
  
  if(!("Cardinal" %in% attr(class(mse),"package"))){
    stop("输入非Cardinal包的数据类型")
  }
  
  if(!is.vector(filename)){
    filename <- as.character(1:length(levels(run(mse))))
  }
  
  if (is.null(layout)) {
    n1 <- ceiling(sqrt(length(filename)))
    n2 <- ceiling(length(filename) / n1)
    layout2 <- c(n2, n1)
  } else {
    layout2 <- layout
    n2 <- layout2[1]
    n1 <- layout2[2]
  }
  
  xmax <- max(coord(mse)$x)
  xmin <- min(coord(mse)$x)
  ymax <- max(coord(mse)$y)
  ymin <- min(coord(mse)$y)
  xlength <- xmax - xmin + 1
  ylength <- ymax - ymin + 1
  
  if(!is.null(figfile)){
    figdata <- png::readPNG(figfile)
  }
  
  if(mapmz){
    if(is.null(mz)){
      realmz <- mz(mse)
    }else{
      realmz <- mz
    }
    
    suppressWarnings(library(foreach))
    suppressWarnings(library(doParallel))
    registerDoParallel(cores=5)
    foreach(j=seq_len(length(realmz))) %dopar% {
    # for (j in seq_len(length(realmz))){
      print(paste0(realmz[[j]],"成像图绘图中"))
      
      if(length(realmz[[j]]) > 1 & !superpose){
        print(realmz[[j]])
        if (is.null(layout)) {
          n1 <- ceiling(sqrt(length(realmz[[j]])))
          n2 <- ceiling(length(realmz[[j]]) / n1)
          layout2 <- c(n2, n1)
        } else {
        }
      }else{
        if (is.null(layout)) {
          n1 <- ceiling(sqrt(length(filename)))
          n2 <- ceiling(length(filename) / n1)
          layout2 <- c(n2, n1)
        } else {
          layout2 <- layout
          n2 <- layout2[1]
          n1 <- layout2[2]
        }
      }
      
      plotfile(savepath = savepath,
               mapname = paste0(mapname, "-",format(realmz[[j]][1], nsmall = 5,trim = T)),
               imagetype = imagetype,
               width = width * n1, height = height * n2,
               family = family,
               units = "px")
      
      if (lightmode) {
        lightmode()
        colorscale <- colorRampPalette(c("white", color))(1000)
        par(bg = "white")
      } else {
        darkmode()
        colorscale <- colorRampPalette(c("black", color))(1000)
      }
      
      showtext::showtext_auto()
      try({
        print(image(mse,
                    mz = realmz[[j]],
                    fun = fun,
                    smooth.image = smooth.image,
                    superpose = superpose,
                    normalize.image = normalize.image,
                    contrast.enhance = contrast.enhance,
                    xlab = xlab,
                    ylab = ylab,
                    colorscale = colorscale,
                    xlim = c(xmin - xlength * 0.05, xmax + xlength * 0.05),
                    ylim = c(ymin - ylength * 0.1, ymax + ylength * 0.1),
                    asp = asp,
                    layout = layout2,
                    zlim = zlim,
                    col = col,
                    ...))
      },silent = F)
      
      if(!is.null(figdata)){
        rasterImage(figdata,
                    xleft = xmin,
                    ybottom = ymax,
                    xright = xmax,
                    ytop = ymin)
      }
      
      plotsave()
      
      showtext::showtext_auto(FALSE)
    }
    
  }else{
    tic <- pixelApply(mse, fun)
    
    plotfile(savepath = savepath,
             mapname = mapname,
             imagetype = imagetype,
             width = width * n1, height = height * n2,
             family = family,
             units = "px")
    
    if (lightmode) {
      lightmode()
      colorscale <- colorRampPalette(c("white", color))(1000)
    } else {
      darkmode()
      colorscale <- colorRampPalette(c("black", color))(1000)
    }
    
    showtext::showtext_auto()
    print(image(mse,
                formula = tic ~ x * y,
                smooth.image = smooth.image,
                superpose = superpose,
                normalize.image = normalize.image,
                contrast.enhance = contrast.enhance,
                xlab = xlab,
                ylab = ylab,
                colorscale = colorscale,
                xlim = c(xmin - xlength * 0.05, xmax + xlength * 0.05),
                ylim = c(ymin - ylength * 0.1, ymax + ylength * 0.1),
                asp = asp,
                subset = subset,
                layout = layout2,
                zlim = zlim,
                ...))
    if(!is.null(figdata)){
      rasterImage(figdata,
                  xleft = xmin,
                  ybottom = ymax,
                  xright = xmax,
                  ytop = ymin)
    }
    
    plotsave()
    showtext::showtext_auto(FALSE)
  }
  
  gc(reset = TRUE)
  setCardinalBPPARAM(SerialParam())
  return("完成绘图")
}

