#' projectpath
#'
#' 根据系统选择项目数据位置
#'
#' @param path 向量，路径
#'
#' @return 返回路径
#' @export
projectpath <- function(path = NULL) {
  # 获取系统
  os <- R.Version()$os
  
  if (grepl("w32", os)) {
    # path1 <- "//192.168.10.177/代谢/lumingos/project/"
    path1 <- "//192.168.10.173/代谢173/lumingos/project-2023/"
  } else if (grepl("linux", os)) {
    # path1 <- "/public/lumingos/project/"
    path1 <- "/data/hstore4/lumingos/project-2023/"
  }
  
  if (!is.null(path)) {
    path1 <- paste0(path1, path)
  }
  
  return(path1)
}

#' projectpath2
#'
#' 项目数据位置
#'
#' @param path 向量，路径
#'
#' @return 以file://返回路径
#' @export
projectpath2 <- function(path = NULL) {
  # 获取系统
  # path1 <- "file://192.168.10.177/代谢/lumingos/project/"
  path1 <- "file://192.168.10.173/代谢173/lumingos/project-2023/"
  if (!is.null(path)) {
    path1 <- paste0(path1, path)
  }
  return(path1)
}
