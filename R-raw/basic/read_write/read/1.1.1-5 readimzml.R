#!/opt/conda/bin/Rscript

#' 读取imzml数据
#'
#' @param filename 文件路径
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param vague 是否虚化
#' @param iteration 虚化次数
#' @param maxPixels 虚化上限
#' @param area 是否绘制选择区域
#' @param areards 选择区域信息路径
#' @param areaname 选择区域名称
#' @param mirrorx 逻辑，x轴镜像对称
#' @param mirrory 逻辑，y轴镜像对称
#' @param trans 逻辑，旋转90
#' @param attach.only 逻辑，数据是否读取到内存
#' @param addy 是否按照y轴添加
#' @param addx 是否按照x轴添加
#' @param filtermz 进行mz筛选
#' @param ... 
#' @param intsenityrange 强度范围
#' @param intsenityratiorange 强度比例范围
#'
#' @export
readimzml <- function(filename,
                      mass.range = NULL,
                      resolution = 0,
                      units = "ppm",
                      vague = F,
                      iteration = 2,
                      maxPixels = 1000000,
                      area = F,
                      areards = NA,
                      areaname = gsub(pattern = "\\.rds$",
                                      replacement = "",
                                      gsub(pattern = "^.*/",replacement = "",x = areards)),
                      mirrorx = rep(F, length(filename)),
                      mirrory = rep(F, length(filename)),
                      trans = rep(F, length(filename)),
                      addy = F,
                      addx = F,
                      filtermz = NULL,
                      intsenityrange = NULL,
                      intsenityratiorange = NULL,
                      attach.only = T,
                      smooth.image = "none",
                      ...) {
  print("~~imzml数据读取开始~~")

  if("package:Cardinal" %in% search()){
    detach("package:Cardinal")
  }
  
  suppressMessages(library("Cardinal"))
  
  if(!is.na(resolution)){
    if(resolution <= 0){
      resolution <- NA
    }
  }
  
  print(paste0("~~",filename[1],"读取中~~"))
  mse <- readMSIData(filename[1],
                     attach.only = attach.only,
                     mass.range = mass.range, resolution = resolution, units = units)
  
  resolution <- resolution(mse)
  
  ibrfilename <- gsub(pattern = "\\.imzML",replacement = ".ibr",x = filename[1])
  if(file.exists(ibrfilename)){
     ibrdata <- readrds(ibrfilename)
     if("pData" %in% names(ibrdata)){
       pData(mse) <- ibrdata$pData
     }
  }
  
  if(is.na(centroided(mse))){
    centroided(mse) <- FALSE
  }
  
  mse$samplename <- gsub(pattern = "-(neg)|-(pos)",replacement = "",x = run(mse))
  
  if (area) {
    print(paste0(areards[1],"读取中"))
    imzmlarea <- readrds(file = areards[1])
    mse <- subsetPixels(mse, imzmlarea)
    run(mse) <- paste(mse$samplename,
                      areaname[1],
                      sep = "-")
    mse$areaname <- areaname[1]
  }else{
    run(mse) <- mse$samplename
  }
  
  if (vague) {
    mse <- imzmlvague(mse = mse,
                      iteration = iteration,
                      maxPixels = maxPixels,
                      smooth.image = smooth.image)
  }
  
  i <- 1
  if (mirrorx[i] & !is.na(mirrorx[i])) {
    coord(mse)$x <- max(coord(mse)$x) + min(coord(mse)$x) - coord(mse)$x
  }
  if (mirrory[i] & !is.na(mirrory[i])) {
    coord(mse)$y <- max(coord(mse)$y) + min(coord(mse)$y) - coord(mse)$y
  }
  
  if (trans[i] & !is.na(trans[i])) {
    coordxy <- coord(mse)$x
    coord(mse)$x <- coord(mse)$y
    coord(mse)$y <- coordxy
  }
  
  if (length(filename) > 1) {
    for (i in 2:length(filename)) {
      print(paste0("~~",filename[i],"读取中~~"))
      mse1 <- readMSIData(filename[i],
                          attach.only = attach.only,
                          mass.range = mass.range, resolution = resolution, units = units)
      
      ibrfilename <- gsub(pattern = "\\.imzML",replacement = ".ibr",x = filename[i])
      if(file.exists(ibrfilename)){
        ibrdata <- readrds(ibrfilename)
        if("pData" %in% names(ibrdata)){
          pData(mse1) <- ibrdata$pData
        }
      }
      
      if(is.na(centroided(mse1))){
        centroided(mse1) <- FALSE
      }
      
      mse1$samplename <- gsub(pattern = "-(neg)|-(pos)",replacement = "",x = run(mse1))
      
      if (area) {
        print(paste0(areards[i],"读取中"))
        imzmlarea <- readrds(file = areards[i])
        mse1 <- subsetPixels(mse1, imzmlarea)
        run(mse1) <- paste(mse1$samplename,
                           areaname[i],
                           sep = "-")
        mse1$areaname <- areaname[i]
      }else{
        run(mse1) <- mse1$samplename
      }
      
      if (vague) {
        mse1 <- imzmlvague(mse = mse1,
                           iteration = iteration,
                           maxPixels = maxPixels)
      }
      
      if (mirrorx[i] & !is.na(mirrorx[i])) {
        coord(mse1)$x <- max(coord(mse1)$x) + min(coord(mse1)$x) - coord(mse1)$x
      }
      if (mirrory[i] & !is.na(mirrory[i])) {
        coord(mse1)$y <- max(coord(mse1)$y) + min(coord(mse1)$y) - coord(mse1)$y
      }
      if (trans[i] & !is.na(trans[i])) {
        coordxy <- coord(mse1)$x
        coord(mse1)$x <- coord(mse1)$y
        coord(mse1)$y <- coordxy
      }
      
      if(addy){
        maxy <- max(coord(mse)$y)
        coord(mse1)$y <- coord(mse1)$y+maxy+20
      }
      
      if(addx){
        maxx <- max(coord(mse)$x)
        coord(mse1)$x <- coord(mse1)$x+maxx+20
      }
      
      # mse3 <<- mse
      # mse4 <<- mse1
      
      mse <- Cardinal::cbind(mse, mse1)
    }
    
  } 
  
  if(is.null(filtermz)){
  }else{
    filtermz <- unique(filtermz)
    filtermz <- as.numeric(filtermz)
    filtermz <- filtermz[order(filtermz)]
    mse <- subsetFeatures(mse,mz = filtermz)
  }
  
  if (!is.null(mass.range)) {
    mse <-  subsetFeatures(mse,mz <= mass.range[2] & mz >= mass.range[1])
  }
  
  if (!is.null(intsenityratiorange)) {
    
    rangedata <- function(x,intsenityratiorange){
      intsenityrange <- intsenityratiorange * max(x)
      x[x > intsenityrange[2]] <- intsenityrange[2]
      x[x < intsenityrange[1]] <- intsenityrange[1]
      return(x)
    }
    
    spectra(mse) <- t(featureApply(mse,rangedata,intsenityratiorange = intsenityratiorange))

  }
  
  if (!is.null(intsenityrange)) {
    spectra(mse)[spectra(mse) > intsenityrange[2]] <- intsenityrange[2]
    spectra(mse)[spectra(mse) < intsenityrange[1]] <- intsenityrange[1]
  }
  
  return(mse)
}

