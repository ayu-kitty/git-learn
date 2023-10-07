#!/opt/conda/bin/Rscript

#' 获取多元统计的模型参数
#'
#' @param data obj 
#'
#' @export
getsummary <- function(rdspath){
  
  data <- readdata(rdspath)
  
  if("mulstatistics" %in% class(data)){
    if(!is.null(data$summarydata)){
      return(data$summarydata)
    }else{
      return(data.frame())
    }
  }else if("mulstatistics-both" %in% class(data)){
    
    lcdata <- getsummary(data$lcresult)
    gcdata <- getsummary(data$gcresult)
    
    data <- rbind(lcdata,gcdata)
    
    return(data)
    
  }else{
    return(data.frame())
  }
}

#' 根据rds获取多元统计的模型参数，并保存
#'
#' @param rdspath rds路径
#' @param savepath 保存路径
#'
#' @export
getsummary_file <- function(rdspath,
                            savepath){

  # 数据提取
  data <- lapply(rdspath,getsummary)

  # 数据处理
  summarydata <- Reduce(rbind,data)
  
  if(dim(summarydata)[1] == 0){
    return()
  }
  
  summarydata <- summarydata[order(summarydata$Q2),]
  summarydata <- summarydata[!duplicated(summarydata[,c("Group","Type")]),]
  summarydata$Type <- factor(summarydata$Type,levels = c("PCA","PLS-DA","OPLS-DA"))
  summarydata <- summarydata[order(summarydata$Type),]
  summarydata <- summarydata[order(summarydata$Group),]
  colnames(summarydata)[colnames(summarydata)=="Type"] <- "Modetype"
  
  # 数据保存
  savetxt(data = summarydata,
          filename = savepath)

}
