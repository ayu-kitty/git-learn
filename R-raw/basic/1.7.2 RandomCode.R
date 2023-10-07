#!/opt/conda/bin/Rscript

#' 随机编码
#'
#' @param n 数值向量，随机编码长度
#'
#' @return 随机字符串
#'
#' @export
RandomCode <- function(n = 20) {
  
  number <- c(as.character(0:9), LETTERS)
  code <- NULL
  
  i <- 1
  while (i <= n) {
    code <- paste0(code, number[floor(runif(1, 1, 36.99))])
    i <- i + 1
  }
  
  return(code)
}
