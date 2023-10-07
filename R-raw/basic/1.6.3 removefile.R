#!/opt/conda/bin/Rscript

#' 删除文件
#'
#' @param path 检索路径
#' @param pattern 文件正则
#'
#' @export
mulremovefile <- function(path = "./",
                          pattern = ".mzML$"){
  filename <- list.files(path = path,
                         pattern = pattern,
                         recursive = T,full.names = T)
  unlink(filename)
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-p","--path",default = "./",
                      help = "检索路径")
  parser$add_argument("-pa","--pattern",default = NULL,
                      help = "文件的正则")
  
  args <- parser$parse_args()
  
  do.call(what = mulremovefile,args = args)
  
}

