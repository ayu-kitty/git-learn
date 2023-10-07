#!/opt/conda/bin/Rscript

#' 获取选区数据及图像
#'
#' @param samplename 样本名称
#' @param samplefrom 样本数据路径
#' @param areaname 区域名
#' @param areafrom 区域数据路径
#' @param mode 正负离子模式
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param asp 图像长宽比
#' @param savemapwd 图片保存路径
#' @param savedatawd 数据保存路径
#' @param imagetype 图片格式
#' @param lightmode 成像模式
#' @param normalize.image 绘图模式
#' @param ... 见[imzmlareaimage()]
#'
#' @export
GetData <- function(samplename,
                    samplefrom = "./sample/final/",
                    areaname,
                    areafrom = "./sample/area/data/",
                    mode,
                    mass.range = NULL,
                    resolution = 5,
                    units = "ppm",
                    asp = 1,
                    savemapwd = "./sample/area/map/",
                    savedatawd = areafrom,
                    imagetype = c("jpg", "pdf"),
                    lightmode = F,
                    normalize.image = "linear",
                    overwrite = T,
                    ...) {
  suppressMessages(library("Cardinal"))
  filename <- paste0(samplefrom, samplename, "-", mode, ".imzML")
  mse <- readMSIData(file = filename, attach.only = F,
                     mass.range = mass.range, resolution = resolution, units = units)
  
  samplearea <- readRDS(file = paste0(areafrom,"/",
                                      samplename,"-", areaname,"-", mode,".rds"))
  
  mz <- data.frame(mz = format(mz(mse), nsmall = 5, trim = T))
  spectradata <- as.data.frame(iData(mse, "intensity")[, samplearea])
  
  titlename <- paste(samplename, areaname,
                     mode, coord(mse[, samplearea])$x, coord(mse[, samplearea])$y,
                     sep = "-")
  names(spectradata) <- titlename
  titlename <- paste(samplename,
                     mode, coord(mse[, samplearea])$x, coord(mse[, samplearea])$y,
                     sep = "-")
  
  meandata <- apply(spectradata, 1, mean)
  meandata <- data.frame(mean = meandata)
  names(meandata) <- paste(samplename, areaname, mode, sep = "-")
  
  spectradata <- cbind(mz, spectradata)
  meandata <- cbind(mz, meandata)
  

  savetxt(data = spectradata,
          filename = paste0(savedatawd,"/",
                            samplename,"-", areaname,"-", mode,"-pixel_level.txt"),
          overwrite = overwrite)
  savetxt(data = meandata,
          filename = paste0(savedatawd,"/",
                            samplename,"-", areaname,"-", mode,"-sample_level.txt"),
          overwrite = overwrite)
  saverds(data = titlename,
          filename = paste0(savedatawd,"/",
                            samplename,"-", areaname,"-", mode,"-name",".rds"),
          overwrite = overwrite)
  
  imzmlareaimage3(filename = filename,
                  areafile = paste0(areafrom,"/",
                                    samplename,"-", areaname,"-", mode,".rds"),
                  savepath = savemapwd,
                  savename = paste0(samplename, "-", areaname, "-", mode),
                  type = imagetype,
                  lightmode = lightmode,
                  normalize.image = normalize.image,
                  asp = asp,
                  ...)
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  parser$add_argument("-a","--asp",default = 1, type= "double",help = "长宽分辨率比")
  parser$add_argument("-s","--samplefrom",default = NULL,help = "样本数据路径")
  
  args <- parser$parse_args()
  
  writeinfo()
  
  createdir(filename = "./sample/area",linkdir = T)
  
  if(is.null(args$samplefrom)){
    if(length(list.files(path = "./sample/final/",pattern = "^qc_data")) > 0){
      args$samplefrom <- "./sample/adjustdata/"
    }else{
      args$samplefrom <- "./sample/final/"
    }
  }
  
  result <- do.call(what = GetAllData,args = args)
  result <- do.call(what = GetClusterData,args = args)
  result <- do.call(what = GetAreaData,args = args)
  result <- do.call(what = GetMulData,args = args)
  
  writeinfo(endtime = T)
  
}