#!/opt/conda/bin/Rscript

#' 批量对目录下imzml根据模块处理
#'
#' @param imzmlpath 数据目录
#' @param moderange 正负离子模式
#' @param imzmlmoudle 调用模块
#' @param ... 
#'
#' @export
imzmlanacall <- function(imzmlpath = "./sample/imzml/",
                         moderange = c("neg", "pos"),
                         imzmlmoudle = NULL,
                         samplename = NULL,
                         ...) {
  # 正负离子循环
  for (mode in moderange) {
    
    # 获取.imzML文件相对路径
    filename <- list.files(path = imzmlpath,
                           pattern = paste0("-", mode, ".imzML$"),
                           full.names = F,
                           recursive = T)
    # 获取玻片名
    slidename <- gsub(pattern = paste0("-", mode, ".imzML"),
                      replacement = "",
                      x = filename)
    
    if(!is.null(samplename)){
      slidename <- slidename[slidename %in% samplename]
    }
    
    # 判断是否找到imzML文件
    if (length(slidename) > 0) {
      
      # 针对imzML文件循环剔除背景处理
      for (i in seq_len(length(slidename))) {
        
        imzmlmoudle(samplename = slidename[i],
                    mode = mode,
                    imzmlpath = imzmlpath,
                    ...)
        gc(reset = TRUE)
        
      }
    } else {
      print(paste0("在",imzmlpath,"目录下未找到", mode, "模式的imzML文件"))
    }
  }
}
