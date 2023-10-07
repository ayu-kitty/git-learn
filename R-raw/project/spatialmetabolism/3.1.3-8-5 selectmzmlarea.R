
#' 获取选区信息
#'
#' @param samplename 样本名称
#' @param areaname 区域名
#' @param mode 正负离子模式
#' @param imzmlpath 数据路径
#' @param saveareawd 区域数据储存位置
#' @param ... 见[imzmlimage()]
#'
#' @export
selectmzmlarea <- function(samplename,
                           areaname,
                           mode = "neg",
                           imzmlpath = "./sample/final/",
                           saveareawd = "./sample/area/data/",
                           ...){
  
  mse <- readdata(filename = paste0(imzmlpath,"/", samplename, "-", mode, ".imzML"))
  imzmlimage(filename = mse,
             imagetype = NA,...)
  samplearea <- selectROI(mse)
  
  saverds(data = samplearea,
          filename = paste0(saveareawd,"/",
                            samplename,"-", areaname,"-", mode,".rds"))
  
}
