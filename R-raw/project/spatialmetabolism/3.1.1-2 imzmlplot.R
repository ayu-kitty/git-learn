#!/opt/conda/bin/Rscript

#' 空代质谱图
#'
#' @param filename 文件路径
#' @param savepath 保存路径
#' @param mapname 保存名称
#' @param imagetype 图片格式
#' @param width 图片宽度
#' @param height 图片长度
#' @param lightmode 成像模式
#' @param superpose 是否分面绘制
#' @param family 字体
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param area 是否绘制选择区域
#' @param areards 选择区域信息路径
#' @param attach.only 逻辑，数据是否读取到内存
#' @param layout 多样本排布方式
#' @param addgroup 逻辑，是否添加组别信息
#' @param mapmz 逻辑值,是否绘制mz的图
#' @param mz mz值，为空时自动提取
#' @param ... 见[Cardinal::plot()]
#'
#' @export
imzmlplot <- function(filename,
                      savepath = "./",
                      mapname = "test",
                      imagetype = c("png","pdf"),
                      mapmz = F,
                      mz = NULL,
                      plusminus = 0.0005,
                      xlim = NULL,
                      width = 8,
                      height = 4,
                      lightmode = T,
                      superpose = F,
                      family = "sans",
                      mass.range = NULL,
                      resolution = 5,
                      units = "ppm",
                      area = F,
                      areards = NA,
                      ylim = NULL,
                      attach.only = F,
                      layout = c(1,1),
                      addgroup = F,
                      fun = mean,
                      intsenityrange = NULL,
                      intsenityratiorange = NULL,
                      ...) {
  print("质谱图绘图开始")

  if("package:Cardinal" %in% search()){
    detach("package:Cardinal")
  }
  suppressMessages(library("Cardinal"))
  
  mse <- readdata(filename = filename,
                  mass.range = mass.range,
                  resolution = resolution,
                  units = units,
                  area = area,
                  areards = areards,
                  intsenityrange = intsenityrange,
                  intsenityratiorange = intsenityratiorange,
                  attach.only = attach.only)
  
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
  
  
  if(mapmz){
    if (length(mz) == 0) {
      realmz <- mz(mse)
    } else {
      realmz <- as.numeric(mz)
    }
    
    for (j in seq_len(length(realmz))) {
      
      mse2 <- mse[mz(mse) < (realmz[j] * (1 + plusminus)),]
      mse2 <- mse2[mz(mse2) > (realmz[j] * (1 - plusminus)),]
      
      if(dim(mse2)[1] == 0){
        mse2 <- mse
        xlim2 <- c(realmz[j] * (1 - plusminus),realmz[j] * (1 + plusminus))
      }else if(dim(mse2)[1] == 1){
        mse3 <- mse
        spectra(mse3)[1, ] <- 0
        msemin <- mse2
        mz(msemin) <- realmz[j] * (1 - plusminus)
        msemax <- mse2
        mz(msemax) <- realmz[j] * (1 + plusminus)
        
        mse2 <- BiocGenerics::rbind(msemin,mse2)
        mse2 <- BiocGenerics::rbind(mse2,msemax)
        spectra(mse2)[c(1, 3), ] <- 0
        xlim2 <- xlim
      }
      
      plotfile(savepath = savepath,
               mapname = paste0(mapname, "-",format(realmz[[j]], nsmall = 5,trim = T)),
               imagetype = imagetype,
               width = width * n1, height = height * n2,
               family = family,
               units = "px")
      
      if (lightmode) { 
        lightmode()
        par(bg = "white")
      } else { darkmode()}
      
      showtext::showtext_auto()
      print(
        if(addgroup){
          pixel.group <<- run(mse)
          Cardinal::plot(mse2,
                         run = run(mse),
                         xlim = xlim,
                         ylim = ylim,
                         superpose = superpose,
                         pixel.groups = pixel.group,
                         layout = c(n2, n1),
                         ...)
        }else{
          Cardinal::plot(mse2,
                         run = run(mse),
                         xlim = xlim,
                         ylim = ylim,
                         superpose = superpose,
                         layout = c(n2, n1),
                         ...)
        }
      )
      
      plotsave()
      showtext::showtext_auto(FALSE)
    }
    
  }else{
    
    plotfile(savepath = savepath,
             mapname = mapname,
             imagetype = imagetype,
             width = width * n1, height = height * n2,
             family = family,
             units = "px")
    
    if (lightmode) { lightmode() } else { darkmode()}
    
    showtext::showtext_auto()
    print(
      if(addgroup){
        pixel.group <<- run(mse)
        Cardinal::plot(mse,
                       superpose = superpose,
                       pixel.groups = pixel.group,
                       xlim = xlim,
                       ylim = ylim,
                       fun = fun,
                       layout = c(n2, n1),...)
      }else{
        Cardinal::plot(mse,
                       superpose = superpose,
                       xlim = xlim,
                       ylim = ylim,
                       fun = fun,
                       layout = c(n2, n1),...)
      }
    )
    
    plotsave()
    
    showtext::showtext_auto(FALSE)
  }
  
  gc(reset = TRUE)
  return("完成绘图")
}

