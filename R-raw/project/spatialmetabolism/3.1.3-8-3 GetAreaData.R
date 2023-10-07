
#' 获取样本选区数据及图像
#'
#' @param samplefrom 样本数据路径
#' @param selectfrom 选区来源路径
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param saveareawd 区域数据路径
#' @param moderange 正负离子模式
#' @param ... 见[GetData()]
#'
#' @export
GetAreaData <- function(mass.range = NULL,
                        resolution = 5,
                        units = "ppm",
                        samplefrom = "./sample/final/",
                        selectfrom = "./sample/msireader",
                        saveareawd = "./sample/area/data/",
                        moderange = c("neg", "pos"),
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
    
    # 判断是否找到imzML文件
    if (length(samplename) > 0) {
      
      # 针对imzML文件循环剔除背景处理
      for (i in seq_len(length(samplename))) {
        areafile <- list.files(path = selectfrom,
                               pattern = paste0("^",samplename[i], "-", ".*-", mode, ".txt$"),
                               full.names = F,
                               recursive = T)
        
        areaname <- gsub(pattern = paste0("-", mode, ".txt"),
                         replacement = "",
                         x = areafile)
        
        areaname <- gsub(pattern = paste0(samplename[i], "-"),
                         replacement = "",
                         x = areaname)
        
        if (length(areaname) > 0) {
          mse <- readMSIData(file = paste0(samplefrom,"/",filename[i]), attach.only = TRUE,
                             mass.range = mass.range, resolution = resolution, units = units)
          
          for (j in seq_len(length(areaname))) {
            print(paste0(samplename[i], "-", areaname[j], "-", mode, ".txt运行中"))
            
            arearange <- read.table(file = paste0(selectfrom,"/",samplename[i],
                                                  "-", areaname[j],
                                                  "-", mode,
                                                  ".txt"))
            
            arearange <- arearange[arearange[, 4] == 1, ]
            
            arearange2 <- paste0(coord(mse)$x,"-", coord(mse)$y) %in% 
              paste0(arearange[, 2],"-", arearange[, 3])
            
            saverds(data = arearange2,
                    filename = paste0(saveareawd,"/",
                                      samplename[i],"-", areaname[j],"-", mode,".rds"))
            
            GetData(samplename = samplename[i],
                    samplefrom = samplefrom,
                    areaname = areaname[j],
                    areafrom = saveareawd,
                    mode = mode,
                    mass.range = mass.range, resolution = resolution, units = units,
                    ...)
          }
        }
      }
    } else {
      print(paste0("在",samplefrom,"目录下未找到", mode, "模式的imzML文件"))
    }
  }
}
