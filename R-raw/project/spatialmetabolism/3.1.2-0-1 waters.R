#!/opt/conda/bin/Rscript

#' 对目录下imzml进行空代数据归一化
#'
#' @param slide 选择玻片名
#' @param imzmlpath 数据路径
#' @param moderange 正负离子模式
#' @param ... 见[getimzmlNormalization()]
#'
#' @export
imzmlwatersdatadeal <- function(slidename = NULL,
                                imzmlpath = "raw",
                                moderange = c("neg", "pos"),
                                savepath = "./sample/imzml/",
                                ...) {
  
  suppressMessages(library("Cardinal"))
  
  # 正负离子循环
  for (mode in moderange) {
    # 数据读取
    if (is.null(slidename)) {
      filename <- list.files(path = paste0(imzmlpath,"/",mode),
                             pattern = "\\.imzML$",
                             full.names = T,
                             recursive = T)
    } else {
      filename <- paste0(imzmlpath,"/",mode,"/",slidename, ".imzML")
      filename <- filename[file.exists(filename)]
    }
    
    if (length(filename) == 0) {
      print(paste0("在", mode, "模式下",imzmlpath,"中未找到imzml文件"))
      next
    }
    # 样本名获取
    samplename <- gsub(pattern = paste0(".imzML"),
                       replacement = "",
                       x = basename(filename))
    
    for ( i in 1:length(samplename)) {
      # 空代数据归一化(见函数getwatersdata2)
      getwatersdata2(filename = filename[i],
                     samplename = samplename[i],
                     mode = mode,
                     savepath = savepath,
                     ...)
      
      # getwatersdata(filename = filename[i],
      #               samplename = samplename[i],
      #               mode = mode,
      #               savepath = savepath,
      #               ...)
      
    }
  }
}

#' 对目录下waters数据归一(对齐版本,已弃用)
#'
#' @export
getwatersdata <- function(filename,
                          samplename,
                          mode,
                          mass.range = NULL,
                          resolution = 20,
                          units = "ppm",
                          tolerance = 0.1, 
                          savepath = "./sample/imzml/",
                          imagetype = "jpg",
                          asp = 1,
                          ...){
  
  suppressMessages(library("Cardinal"))
  # setCardinalBPPARAM(MulticoreParam(workers = 3))
  
  if(!is.na(tolerance)){
    if(tolerance <= 0){
      tolerance <- 0.1
    }
  }else{
    tolerance <- 0.1
  }
  
  mse <- readdata(filename = filename,
                  mass.range = mass.range,
                  resolution = resolution,
                  units = units,
                  ...)
  
  print(mse)
  print(paste0("~使用分辨率为:",resolution(mse),names(resolution(mse))))
  print(paste0("~使用峰校正参数为:",tolerance,"mz"))
  
  if(mode == "neg"){
    mz <- 554.26202254
  }else{
    # mz <- c(104.107539075,556.27657454)
    mz <- 556.27657454
  }
  
  temppath <- NULL
  
  source(file = packagepath(path = "script/spatial/mzAlign.R"))
  
  if(length(pixels(mse)) > 50000){
    
    temppath <- "temp"
    
    maxnum <- 30000
    
    num <- round(length(pixels(mse))/maxnum)
    filename2 <- NULL
    
    i <- 1
    
    for ( i in 1:num) {
      
      base::print(i)
      
      startnum <- (i-1)*maxnum+1
      if(i == num){
        endnum <- length(pixels(mse))
      }else{
        endnum <- i*maxnum
      }
      
      mse3 <- subsetPixels(mse,startnum:endnum)
      
      mse2 <- mzAlign2(mse3, ref = mz,tolerance = tolerance, units = "mz") %>% process()
      
      mse2 <- as(mse2, "MSProcessedImagingExperiment")
      # resolutiondata <- c(resolution)
      # names(resolutiondata) <- units
      # resolution(mse2) <- resolutiondata
      # mse2 <- as(mse2, "MSContinuousImagingExperiment")
      
      # print(mse2)
      
      # mse2 <- Cardinal::pull(mse2, as.matrix = TRUE)
      # data <- as.matrix(spectra(mse2))
      # data[data < 1] <- 0
      # spectra(mse2) <- data
      
      filename <- paste0(temppath,"/",samplename,"-",mode,"-",i,".imzML")
      filename2 <- c(filename2,filename)
      
      saveimzml(data = mse2,
                filename = filename)
      gc(reset = TRUE)
    }
    
    mse2 <- readdata(filename = filename2,
                     mass.range = mass.range,
                     resolution = resolution,
                     units = units,
                     ...)
    
  }else{
    mse2 <- mzAlign2(mse, ref = mz,tolerance = 0.1, units = "mz") %>% process()
    
    mse2 <- as(mse2, "MSProcessedImagingExperiment")
    class(mse2) <- class(mse)
    # resolutiondata <- c(resolution)
    # names(resolutiondata) <- units
    # resolution(mse2) <- resolutiondata
    # mse2 <- as(mse2, "MSContinuousImagingExperiment")
    
    # mse2 <- pull(mse2, as.matrix = TRUE)
    # data <- as.matrix(spectra(mse2))
    # data[data < 1] <- 0
    # spectra(mse2) <- data
  }
  
  filename <- paste0(savepath,"/",samplename,"-",mode,".imzML")
  
  saveimzml(data = mse2,
            filename = filename,
            intensity.type = "64-bit integer")
  
  imzmlplot(filename = mse2,
            mapmz = T,
            mz = mz,
            savepath = savepath,
            mapname = paste0(samplename,"-",mode),
            imagetype = imagetype)
  
  imzmlimage(filename = mse2,
             savepath = savepath,
             mapname = paste0(samplename,"-",mode),
             imagetype = imagetype,
             asp = asp)
  
  if(!is.null(temppath)){
    unlink(temppath,recursive = T,force = T)
  }
  
  # setCardinalBPPARAM(SerialParam())
}


#' 对目录下imzml进行空代数据归一化-waters
#'
#' @param filename 选择玻片名
#' @param samplename 选择样本名
#' @param mode 离子模式
#' @param mass.range mz范围
#' @param resolution 分辨率(此处无用传参)
#' @param units 分辨率单位(此处无用传参)
#' @param tolerance waters峰校正参数,默认20(此处无用传参)
#' @param savepath 保存路径
#' @param imagetype 图片格式(内标mz质谱图)
#' @param asp 图像长宽比
#' 
#' @export
getwatersdata2 <- function(filename,
                           samplename,
                           mode,
                           mass.range = NULL,
                           resolution = 20,
                           units = "ppm",
                           tolerance = 0.1,
                           savepath = "./sample/imzml/",
                           imagetype = "jpg",
                           asp = 1,
                           ...){
  
  suppressMessages(library("Cardinal"))
  # setCardinalBPPARAM(MulticoreParam(workers = 3))
  
  if(!is.na(tolerance)){
    if(tolerance <= 0){
      tolerance <- NA
    }
  }
  # 读数据
  mse <- readdata(filename = filename,
                  mass.range = mass.range,
                  resolution = resolution,
                  units = units,
                  ...)
  
  print(mse)
  
  filename <- paste0(savepath,"/",samplename,"-",mode,".imzML")
  # 转存到imzml文件夹
  saveimzml(data = mse,
            filename = filename)
  #画成像图
  imzmlimage(filename = mse,
             savepath = savepath,
             mapname = paste0(samplename,"-image-",mode),
             imagetype = imagetype,
             asp = asp)
  # 找到内标物质(可以考虑参数释放)
  if(mode == "neg"){
    mz <- 554.26202254
  }else{
    # mz <- c(104.107539075,556.27657454)
    mz <- 556.27657454
  }
  # 内标物质质谱图
  imzmlplot(filename = mse,
            mapmz = T,
            mz = mz,
            savepath = savepath,
            mapname = paste0(samplename,"-plot-",mode),
            imagetype = imagetype)
  # 获取基础的imzml info 信息(函数在1.3中)
  getimzmlbasicinfo(mse = mse,
                    mode = mode,
                    waters = T,
                    filename = paste0(savepath,"/",samplename, "-", mode,".xlsx"))
  getimzmlbasicinfo(mse = mse[,1:10],
                    mode = mode,
                    waters = T,
                    filename = paste0(savepath,"/",samplename, "-", mode,"-bg.xlsx"))
  
  # setCardinalBPPARAM(SerialParam())
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-i","--imzmlpath",default = "./raw/", help = "imzml原始文件路径,默认./raw/")
  parser$add_argument("-s","--savepath",default = "./sample/imzml/", help = "imzml数据保存路径,默认./sample/imzml/")
  parser$add_argument("-m","--massrange",default = NULL,type= "double",help = "mz范围,默认NULL",
                      nargs = 2,dest = "mass.range")
  parser$add_argument("-re","--resolution",default = 20, type= "double",help = "分辨率,默认20")
  parser$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  parser$add_argument("-tr","--tolerance",default = 0.1, type= "double",help = "峰校正参数,默认0.1")
  parser$add_argument("-a","--asp",default = 1, type= "double",help = "长宽比")
  args <- parser$parse_args()
  
  writeinfo()
  createdir(filename = args$savepath,linkdir = T)
  
  mulargs <- do.call(what = imzmlwatersdatadeal,args = args)
  
  writeinfo(endtime = T)
  
}
