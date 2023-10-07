
#' 返回对应向量,同match
#'
#' @param a 两列dataframe
#' @param b 向量
#'
#' @return 将a中第一列与b对应的第二列向量
#'
#' @export
whto <- function(a, b) {
  
  if(any(duplicated(a[,1]))){
    warning("对应矩阵有重复",immediate. = T)
  }
  
  if(!all(b %in% a[,1])){
    stop("对应矩阵有缺失")
  }
  
  group <- a[match(b,a[,1]),2]

  return(group)
}