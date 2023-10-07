#!/opt/conda/bin/Rscript

#' 获取多元统计的loading数据
#'
#' @param data obj
#'
#' @export
getloading <- function(data) {
  
  if(!is.null(data$statistics)){
    statistics <- data$statistics
    mode <- statistics@typeC
    if (mode == "PCA") {
      score <- data.frame(num = row.names(statistics@loadingMN),
                          statistics@loadingMN[, 1:3],
                          check.names = F, stringsAsFactors = F)
      
      x <- paste0("PC1(", statistics@modelDF[1, 1] * 100, "%)")
      y <- paste0("PC2(", statistics@modelDF[2, 1] * 100, "%)")
      z <- paste0("PC3(", statistics@modelDF[3, 1] * 100, "%)")
    }else if(mode == "PLS-DA"){
      score <- data.frame(num = row.names(statistics@loadingMN),
                          statistics@loadingMN[, 1:3],
                          VIP = statistics@vipVn,
                          check.names = F, stringsAsFactors = F)
      
      x <- paste0("PC1(", statistics@modelDF[1, 1] * 100, "%)")
      y <- paste0("PC2(", statistics@modelDF[2, 1] * 100, "%)")
      z <- paste0("PC3(", statistics@modelDF[3, 1] * 100, "%)")
    }else {
      score <- data.frame(num = row.names(statistics@loadingMN),
                          statistics@loadingMN[, 1, drop = F],
                          statistics@orthoLoadingMN[, 1:2, drop = F],
                          VIP = statistics@vipVn,
                          check.names = F, stringsAsFactors = F)
      x <- paste0("PC1(", statistics@modelDF[1, 1] * 100, "%)")
      y <- paste0("PCo1(", statistics@modelDF[2, 1] * 100, "%)")
      z <- paste0("PCo1(", statistics@modelDF[3, 1] * 100, "%)")
    }
    
    colnames(score)[1:4] <- c("num",x,y,z)
    
  }else{
    score <- NULL
  }
  
  return(score)
}

#' 根据rds获取多元统计的loading数据，并保存
#'
#' @param rdspath rds路径
#' @param infofile info信息文件路径
#' @param savepath 保存路径
#'
#' @export
getloading_file <- function(rdspath,
                            infofile = NULL,
                            savepath = NULL){

  # 数据提取
  data <- readdata(rdspath)
  info <- readdata(infofile,row.names = 1)

  # 数据处理
  # mode <- data$statistics@typeC
  score <- getloading(data)
  
  if(is.null(score)){
    return()
  }
  
  if(is.null(info)){
    colnames(score)[1] <- "ID"
  }else{
    score <- merge(score,info,by.x="num",by.y=0,sort=F)
    if(!("ID" %in% colnames(score))){
      colnames(score)[1] <- "ID"
    }
  }
  
  if(is.null(savepath)){
    return(score)
  }else{
    # 数据保存
    savetxt(data = score,
            filename = savepath)
  }

}
