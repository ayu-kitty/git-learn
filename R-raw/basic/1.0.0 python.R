#!/opt/conda/bin/Rscript

#' 加载python
#'
#' @export
loadpython <- function() {
  suppressMessages(library("reticulate"))
  use_condaenv(condaenv = "/opt/conda")
  py_available(initialize = TRUE)
}

#' 运行python脚本
#'
#' @param path python路径
#'
#' @export
runpython <- function(path) {
  loadpython()
  source_python(file = path)
}
