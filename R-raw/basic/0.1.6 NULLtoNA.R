#' NULLtoNA
#'
#' 将NULL转换为NA或其他字符
#'
#' @param data 向量或矩阵
#' @param nastop 逻辑，为空时报错
#' @param fill 向量，填充字符，默认为NA
#'
#' @export
NULLtoNA <- function(data,
                     nastop = F,
                     fill = NA) {
  if (is.null(data)) {
    if (nastop) {
      stop("关键索引为空，停止运行")
    }
    return(fill)
  } else if (length(data) == 0) {
    if (nastop) {
      stop("关键索引为空，停止运行")
    }
    if (is.vector(data)) {
      return(fill)
    } else if (is.data.frame(data) | is.matrix(data)) {
      if (nrow(data) == 0) {
        return(fill)
      } else {
        return(rep(fill, nrow(data)))
      }
    } else {
      warning("数据格式没查询到，返回原数据", immediate. = T)
      return(data)
    }
  } else {
    return(data)
  }
}
