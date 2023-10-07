
#' 质谱成像图
#'
#' @param samplename 样本名
#' @param imagename 图像名称
#' @param area 是否绘制选择区域
#' @param areaname 区域名
#' @param areapath 区域信息文件路径
#' @param mzmlpath 数据路径
#' @param mode 正负离子模式
#' @param imagetype 图片格式
#' @param savepath 保存路径
#' @param mz mz值，为空时自动提取
#' @param mapmoudle 绘图模块
#' @param ananame 保存文件名
#' @param mapname 图片名
#' @param ... 见mapmoudle调用函数
#'
#' @export
imzmlimagemap <- function(samplename = NULL,
                          imagename = samplename[1],
                          mode = "neg",
                          ananame = paste0(imagename,"/",mode),
                          area = F,
                          areapath = "./area/data/",
                          areaname = "ALL",
                          imzmlpath = "./sample/final/",
                          savepath = "./sample/map/Intensity/",
                          imagetype = c("jpg", "pdf"),
                          mapmoudle = imzmlimage,
                          mapname = "Intensity-",
                          ...) {
  print(paste0(samplename, "-", mode, "运行中"))
  filename <- paste0(imzmlpath, samplename, "-", mode, ".imzML")
  areards <- paste0(areapath,samplename,"-",areaname,"-", mode, ".rds")
  
  mapmoudle(filename = filename,
            area = area,
            areards = areards,
            savepath =   ifelse(is.null(ananame),
                                savepath,
                                paste0(savepath,"/",ananame,"/")),
            mapname = paste0(mapname, imagename, "-", mode),
            imagetype = imagetype,
            ...)
  
}

#' 批量质谱成像图
#'
#' @param area 是否绘制选择区域
#' @param areaname 区域名
#' @param areapath 区域信息文件路径
#' @param imzmlpath 数据路径
#' @param mulname 多样本名称
#' @param moderange 正负离子模式
#' @param singleanalyst 逻辑，是否单个样本分析
#' @param allanalyst 逻辑，是否多个样本分析
#' @param ... 见[imagemap()]
#'
#' @export
Mulimagemap <- function(imzmlpath = "./sample/final/" ,
                        area = F,
                        areaname = "ALL",
                        areapath = "./area/data/",
                        mulname = "all",
                        moderange = c("neg", "pos"),
                        delsample = c("bg_data","qc_data"),
                        singleanalyst = T,
                        allanalyst = T,
                        ...) {
  
  if(is.null(imzmlpath)){
    if(length(list.files(path = "./sample/final/",pattern = "^qc_data")) > 0){
      imzmlpath <- "./sample/adjustdata/"
    }else{
      imzmlpath <- "./sample/final/"
    }
  }
  
  for (mode in moderange) {
    
    # 获取.imzML文件绝对路径
    filename <- list.files(path = imzmlpath,
                           pattern = paste0("-", mode, ".imzML$"),
                           full.names = F,
                           recursive = T)
    if(length(filename) == 0){
      print(paste0("在",imzmlpath,"目录下未找到", mode, "模式的imzML文件"))
      next
    }
    
    # 获取样品名
    samplename <- gsub(pattern = paste0("-", mode, ".imzML"),
                       replacement = "",
                       x = filename)
    if(!is.null(delsample)){
      for(delsample2 in delsample) {
        samplename <- samplename[!grepl(pattern = paste0("^",delsample2,".*"),x = samplename)]
      }
    }
    imagename <- samplename
    
    areaname2 <- NULL
    samplename2 <- NULL
    if(area){
      for ( i in seq_len(length(samplename))) {
        areards <- paste0(areapath,samplename[i],"-",areaname,"-", mode, ".rds")
        areaname3 <- areaname[file.exists(areards)]
        samplename3 <- rep(samplename[i],length(areaname3))
        areaname2 <- c(areaname2,areaname3)
        samplename2 <- c(samplename2,samplename3)
      }
      samplename <- samplename2
      imagename <- paste0(samplename,"-",areaname2)
    }
    
    # 判断是否找到imzML文件
    if (length(samplename) > 0) {
      if(singleanalyst){
        # 针对imzML文件循环剔除背景处理
        for (i in seq_len(length(samplename))) {
          
          imzmlimagemap(samplename = samplename[i],
                        imagename = imagename[i],
                        imzmlpath = imzmlpath,
                        mode = mode,
                        area = area,
                        areaname = areaname2[i],
                        areapath =  areapath,
                        ...)
          gc(reset = TRUE)
          
        }
      }
      
      if(allanalyst){
        if (length(samplename) > 1) {
          
          imzmlimagemap(samplename = samplename,
                        imagename = mulname,
                        imzmlpath = imzmlpath,
                        mode = mode,
                        area = area,
                        areaname = areaname2,
                        areapath =  areapath,
                        ...)
          gc(reset = TRUE)
          
        }
      }else{
        print(paste0("在",imzmlpath,"目录下仅找到一个", mode, "模式的imzML文件"))
      }
      
    } else {
      print(paste0("在",imzmlpath,"目录下未找到", mode, "模式的imzML文件"))
    }
  }
}

#' 批量绘制累加质谱成像图
#'
#' @param mapmoudle 绘图模块
#' @param ananame 保存文件名
#' @param ... 见[Mulimagemap()]
#'
#' @export
getsumimzmlimage <- function(mapmoudle = imzmlimage,
                             ananame = NULL,
                             ...){
  Mulimagemap(mapmoudle = mapmoudle,
              ananame = ananame,
              ...)
}


#' 批量绘制表达物质数的质谱成像图
#'
#' @param mapmoudle 绘图模块
#' @param ananame 保存文件名
#' @param mapname 图片名
#' @param fun 调用函数
#' @param smooth.image 平滑模式
#' @param savepath 保存路径
#' @param ... 见[Mulimagemap()]
#'
#' @export
getnumberimzmlimage <- function(mapmoudle = imzmlimage,
                                ananame = NULL,
                                mapname = "Number-",
                                fun = calzeronumber,
                                smooth.image = "none",
                                savepath = "./sample/map/",
                                ...){
  Mulimagemap(mapmoudle = mapmoudle,
              ananame = ananame,
              mapname = mapname,
              fun = fun,
              smooth.image = smooth.image,
              savepath = savepath,
              ...)
}

