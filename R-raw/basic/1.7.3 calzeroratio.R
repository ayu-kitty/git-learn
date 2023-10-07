#!/opt/conda/bin/Rscript

#' 数值比例计算
#'
#' @param x 数值向量
#' @param intensity 筛选强度
#'
#' @export
calzeroratio <- function(x,intensity = 0){
  
  sum(x > intensity)/length(x)
  
}
