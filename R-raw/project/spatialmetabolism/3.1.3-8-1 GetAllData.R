
#' 获取整个样本数据及图像
#'
#' @param samplefrom 样本数据路径
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param saveareawd 区域数据路径
#' @param moderange 正负离子模式
#' @param ... 见[GetData()]
#'
#' @export
GetAllData <- function(mass.range = NULL,
                       resolution = 5,
                       units = "ppm",
                       samplefrom = "./sample/final/",
                       saveareawd = "./sample/area/data/",
                       moderange = c("neg", "pos"),
                       delsample = c("bg_data","qc_data"),
                       ...) {
  suppressMessages(library("Cardinal"))
  
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
    
    if(!is.null(delsample)){
      for(delsample2 in delsample) {
        samplename <- samplename[!grepl(pattern = paste0("^",delsample2,".*"),x = samplename)]
      }
    }
    
    # 判断是否找到imzML文件
    if (length(samplename) > 0) {
      
      # 针对imzML文件循环剔除背景处理
      for (i in seq_len(length(samplename))) {
        
        mse <- readMSIData(file = paste0(samplefrom,"/",filename[i]), 
                           attach.only = TRUE,
                           mass.range = mass.range, 
                           resolution = resolution, 
                           units = units)
        
        arearange2 <- rep(T, length(mse))
        
        saverds(data = arearange2,
                filename = paste0(saveareawd,"/",
                                  samplename[i],"-", "ALL", "-", mode,".rds"))
        
        GetData(samplename = samplename[i],
                samplefrom = samplefrom,
                areaname = "ALL",
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