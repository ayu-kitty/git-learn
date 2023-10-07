#!/opt/conda/bin/Rscript

#' 获取多元统计的splot数据
#'
#' @param data obj
#'
#' @export
getsplot <- function(data) {
  
  if(!is.null(data$statistics)){
    statistics <- data$statistics
    mode <- statistics@typeC
    if (mode != "PCA") {
      tCompMN <- statistics@scoreMN[, 1, drop = F]
      cxtCompMN <- cor(statistics@suppLs[["xModelMN"]],
                       tCompMN,
                       use = "pairwise.complete.obs")
      score <- data.frame(num = row.names(statistics@loadingMN),
                          p1 = statistics@loadingMN[, 1],
                          p2 = cxtCompMN[, 1],
                          vip = statistics@vipVn,
                          check.names = F, stringsAsFactors = F)
      
      x <- paste0("PC1(", statistics@modelDF[1, 1] * 100, "%)")
      y <- "PC1(corr)"
      colnames(score) <- c("num",x, y, "VIP")
      
    }else{
      score <- NULL
    }
    
  }else{
    score <- NULL
  }
  
  return(score)
}

#' 根据rds获取多元统计的splot数据，并保存
#'
#' @param rdspath rds路径
#' @param infofile info信息文件路径
#' @param savepath 保存路径
#'
#' @export
getsplot_file <- function(rdspath,
                          infofile = NULL,
                          savepath = NULL){

  # 数据提取
  data <- readdata(rdspath)
  info <- readdata(infofile,row.names = 1)

  # 数据处理
  # mode <- data$statistics@typeC
  score <- getsplot(data)
  
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
