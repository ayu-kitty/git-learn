#!/opt/conda/bin/Rscript

#' 加载当前R脚本路径的脚本
#' 
#' @param filename 加载文件名
#'
#' @export
source_script <- function(filename,
                          ...){
  data_try <- try({
    filepath <- paste0(dirname(whereami::thisfile()),"/")
  })
  
  if ("try-error" %in% class(data_try)) {
    filepath <- "./"
  }
  
  filename <- paste0(filepath,filename)
  print("加载",filename,"文件")
  source(filename,...)
}
