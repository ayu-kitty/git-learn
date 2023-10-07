
#' 获取多区域并集与差集
#'
#' @param samplefrom 样本数据路径
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param saveareawd 区域数据路径
#' @param selectpath 选区来源路径
#' @param infopath 项目登记单路径
#' @param moderange 正负离子模式
#' @param ... 见[GetData()]
#' 
#' @export
GetMulData <- function(mass.range = NULL,
                       resolution = 5,
                       units = "ppm",
                       samplefrom = "./sample/final/",
                       saveareawd = "./sample/area/data/",
                       moderange = c("neg", "pos"),
                       selectpath = "./sample/area/data/",
                       infopath = "项目登记单.xlsx",
                       ...) {
  rawareadata <- readdata(filename = infopath, sheet = "多区域选区")
  
  if(is.null(rawareadata)){
    print("项目登记单中无多区域选区信息")
    return()
  }else if(dim(rawareadata)[1] == 0) {
    print("项目登记单中无多区域选区信息")
    return()
  }
  
  # 正负离子循环
  for (mode in moderange) {
    areadata <- rawareadata[rawareadata$模式 == mode | rawareadata$模式 == "both", ]
    print(areadata)
    
    # 判断是否找到imzML文件
    if (dim(areadata)[1] > 0) {
      
      # 针对imzML文件循环剔除背景处理
      for (i in 1:dim(areadata)[1]) {
        print(i)
        
        samplename <- areadata$样品[i]
        needarea <- unlist(strsplit(x = areadata$交集[i], split = "\\+"))
        needdata <- IntersectionArea(samplename = samplename,
                                     needarea = needarea,
                                     mode = mode,
                                     selectpath = selectpath)
        
        if (!is.na(areadata$差集[i])) {
          cutarea <- unlist(strsplit(x = areadata$差集[i], split = "\\+"))
          cutdata <- IntersectionArea(samplename = samplename,
                                      needarea = cutarea,
                                      mode = mode,
                                      selectpath = selectpath)
          
          needdata[cutdata] <- F
        }
        
        areaname <- areadata$选区[i]
        
        saverds(data = needdata,
                filename = paste0(saveareawd,"/",
                                  samplename,"-", areaname,"-", mode,".rds"))
        
        GetData(samplename = samplename,
                samplefrom = samplefrom,
                areaname = areaname,
                areafrom = saveareawd,
                mode = mode,
                mass.range = mass.range, resolution = resolution, units = units,
                ...)
        
      }
    } else {
      print(paste0("在项目登记单下未找到", mode, "模式的多区域选区"))
    }
  }
}


#' IntersectionArea
#'
#' 生成imzmlreader的并集文件
#'
#' @param samplename 样本名
#' @param needarea 区域名
#' @param mode 正负离子模式
#' @param selectpath 选区来源路径
#'
#' @export
IntersectionArea <- function(samplename,
                             needarea,
                             mode,
                             selectpath = "./sample/area/data/") {
  arearange2 <- NULL
  
  for (j in seq_len(length(needarea))) {
    print(paste0(samplename, "-", needarea[j], "-", mode, ".rds数据提取中"))
    arearange <- readdata(filename = paste0(selectpath, "/",
                                            samplename,"-", needarea[j],"-", mode,".rds"))
    if(is.null(arearange2)){
      arearange2 <- arearange
    }else{
      arearange2 <- arearange2|arearange
    }
    
  }
  
  return(arearange2)
}
