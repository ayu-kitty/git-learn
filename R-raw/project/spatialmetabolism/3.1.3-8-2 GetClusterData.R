
#' 根据cluster获取样本区域数据及图像
#'
#' @param samplefrom 样本数据路径
#' @param clusterfrom cluster数据路径
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param saveareawd 区域数据路径
#' @param moderange 正负离子模式
#' @param ... 见[GetData()]
#'
#' @export
GetClusterData <- function(mass.range = NULL,
                           resolution = 5,
                           units = "ppm",
                           samplefrom = "./sample/final/",
                           clusterfrom = "./sample/cluster/sscc/",
                           saveareawd = "./sample/area/data/",
                           moderange = c("neg", "pos"),
                           ...) {
  suppressMessages(library("Cardinal"))
  
  rawareadata <- readdata(filename = "项目登记单.xlsx", sheet = "聚类选区")
  
  if(is.null(rawareadata)){
    print("项目登记单中无多区域选区信息")
    return()
  }else if(dim(rawareadata)[1] == 0) {
    print("项目登记单中无多区域选区信息")
    return()
  }
  
  # 正负离子循环
  for (mode in moderange) {
    
    # 获取.imzML文件相对路径
    filename <- list.files(path = samplefrom,
                           pattern = paste0("-", mode, ".imzML$"),
                           full.names = F,
                           recursive = T)
    
    # 获取玻片名
    samplename <- gsub(pattern = paste0("-", mode, ".imzML"),
                       replacement = "",
                       x = filename)
    
    areadata <- rawareadata[rawareadata$模式 == mode, ]
    areadata <- areadata[areadata$样品 %in% samplename, ]
    print(areadata)
    
    # 判断是否找到imzML文件
    if (dim(areadata)[1] > 0) {
      
      # 针对imzML文件循环剔除背景处理
      for (i in 1:dim(areadata)[1]) {
        
        print(areadata[i,])
        if (areadata[i, "聚类来源"] == "all") {
          ssc <- readdata(filename = paste0(clusterfrom,"/all-", mode,"-data.rds"))
          ssc <- ssc$sscc
          ssccluster <- data.frame(ssc@elementMetadata@coord,class = ssc$class[[1]],sample = run(ssc),select = T)
          ssccluster <- ssccluster[ssccluster$sample == paste0(areadata[i, "样品"],"-", mode),]
        } else if (areadata[i, "聚类来源"] == "single") {
          ssc <- readdata(filename = paste0(clusterfrom,"/",areadata[i, "样品"],"-", mode, "-data.rds"))
          ssc <- ssc$sscc
          ssccluster <- data.frame(ssc@elementMetadata@coord,class = ssc$class[[1]],sample = run(ssc),select = T)
        }
        
        # arearange <- areadata[i, "聚类"]
        arearange <- as.character(areadata[i, "聚类"])
        arearange <- as.numeric(strsplit(x = arearange, split = "\\+")[[1]])
        ssccluster <- ssccluster[ssccluster$class %in% arearange,]
        
        mse <- readMSIData(file = paste0(samplefrom,areadata[i, "样品"],"-",mode,".imzML"),
                           attach.only = TRUE,
                           mass.range = mass.range, resolution = resolution, units = units)
        msearea <- data.frame(coord(mse),select = F)
        msearea[paste(msearea$x,msearea$y) %in% paste(ssccluster$x,ssccluster$y),"select"] <- T
        
        arearange2 <- msearea$select
        
        saverds(data = arearange2,
                filename = paste0(saveareawd,"/",
                                  areadata[i, "样品"],"-", areadata[i, "选区"],"-", mode,".rds"))
        
        GetData(samplename = areadata[i, "样品"],
                samplefrom = samplefrom,
                areaname = areadata[i, "选区"],
                areafrom = saveareawd,
                mode = mode,
                mass.range = mass.range, resolution = resolution, units = units,
                ...)
        
      }
    } else {
      print(paste0("在",samplefrom,"目录下未找到", mode, "模式的imzML文件"))
    }
  }
}
