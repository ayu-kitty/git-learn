#!/opt/conda/bin/Rscript

#' 获取预处理结果数据
#' 
#' @param data 数据
#' @param class 样本分组
#' @param trans 逻辑值，是否转置
#' @param filename 保存文件名
#'
#' @export
readpredata <- function(filename = "data_normalized.csv",
                        colname){
  
  predata <- readdata(filename = filename,
                      row.names = 1)
  if(row.names(predata)[1] == "Label"){
    
    predata <- predata[-1,]
    predata[] <- apply(t(predata), 1, as.numeric)
    
  }else{
    predata <- predata[,-1]
    predata <- data.frame(t(predata),
                          check.names = F,
                          stringsAsFactors = F)
  }
  
  predata <- predata[,colname,drop = F]
  
  return(predata)
}