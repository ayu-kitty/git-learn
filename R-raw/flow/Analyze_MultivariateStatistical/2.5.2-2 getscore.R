#!/opt/conda/bin/Rscript

#' 获取多元统计的score数据
#'
#' @param data obj
#'
#' @export
getscore <- function(data) {
  
  if(!is.null(data$statistics)){
    statistics <- data$statistics
    mode <- statistics@typeC
    if (mode == "PCA" | mode == "PLS-DA") {
      score <- data.frame(statistics@scoreMN[, 1:3],
                          check.names = F, stringsAsFactors = F)
      
      x <- paste0("PC1(", statistics@modelDF[1, 1] * 100, "%)")
      y <- paste0("PC2(", statistics@modelDF[2, 1] * 100, "%)")
      z <- paste0("PC3(", statistics@modelDF[3, 1] * 100, "%)")
    } else {
      score <- data.frame(cbind(statistics@scoreMN[, 1, drop = F],
                                statistics@orthoScoreMN[, 1:2, drop = F]),
                          check.names = F, stringsAsFactors = F)
      
      x <- paste0("PC1(", statistics@modelDF[1, 1] * 100, "%)")
      y <- paste0("PCo1(", statistics@modelDF[2, 1] * 100, "%)")
      z <- paste0("PCo2(", statistics@modelDF[3, 1] * 100, "%)")
    }
    
    colnames(score) <- c(x,y,z)
    score <- merge(score,data$class,by = 0,sort=F)
    colnames(score)[1] <- "Sample"
  }else{
    score <- NULL
  }
  
  return(score)
}

#' 根据rds获取多元统计的score数据，并保存
#'
#' @param rdspath rds路径
#' @param savepath 保存路径
#'
#' @export
getscore_file <- function(rdspath,
                          savepath = NULL){

  # 数据提取
  data <- readdata(rdspath)

  # 数据处理
  # mode <- data$statistics@typeC
  score <- getscore(data)
  
  if(is.null(score)){
    return()
  }

  if(is.null(savepath)){
    return(score)
  }else{
    # 数据保存
    savetxt(data = score,
            filename = savepath)
  }

}
