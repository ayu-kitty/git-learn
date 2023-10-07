#!/opt/conda/bin/Rscript

#' 异常值填NA
#'
#' @param x 数值向量
#' @param time.iqr 异常值范围
#'
#' @export
fun.outlier <- function(x,
                        time.iqr = 3) {
  # if(sum(x!=0) > 20){
  #   x <- x[x!=0]
  #   if(length(x) == 0){
  #     return(0)
  #   }
  # }
  outlier.low <- quantile(x,probs=c(0.25))-IQR(x)*time.iqr
  outlier.high <- quantile(x,probs=c(0.75))+IQR(x)*time.iqr
  x[which(x>outlier.high | x<outlier.low)] <- NA
  return(x)
}

#' 异常值去除后的最大值
#'
#' @param x 数值向量
#' @param time.iqr 异常值范围
#'
#' @export
outliermax <- function(x,
                       time.iqr = 3,
                       ...){
  x <- fun.outlier(x = x,
                   time.iqr = time.iqr)
  max(x,na.rm = T,...)
}

#' 异常值去除后的均值
#'
#' @param x 数值向量
#' @param time.iqr 异常值范围
#'
#' @export
outliermean <- function(x,
                        time.iqr = 3,
                        ...){
  x <- fun.outlier(x = x,
                   time.iqr = time.iqr)
  mean(x,na.rm = T,...)
}


#' 用去除异常值后的最大值填充异常值
#'
#' @param x 数值向量
#' @param time.iqr 异常值范围
#'
#' @export
outliertomax <- function(x,
                         time.iqr = 3) {
  data <- x
  if(sum(x!=0) > 20){
    x <- x[x!=0]
  }else{
    return(data)
  }
  outlier.high <- quantile(x,probs=c(0.75))+IQR(x)*time.iqr
  maxdata <- max(x[which(x<outlier.high)])
  data[data > outlier.high] <- maxdata
  return(data)
}
