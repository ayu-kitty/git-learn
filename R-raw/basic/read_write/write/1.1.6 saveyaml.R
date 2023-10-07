#!/opt/conda/bin/Rscript

#' 保存yaml文件
#'
#' @param data list数据
#' @param filename 文件名
#'
#' @export
saveyaml <- function(data,
                     filename,
                     overwrite = T,
                     ...) {
  
  suppressMessages(library("yaml"))

  if(is.null(data)) {return()}
  if(!is.list(data)) {return()}
  
  print(paste0("在", filename, "中保存"))
  
  dirfilename <- createdir(filename = dirname(filename))
  
  if(!file.exists(filename) | overwrite){
    write_yaml(x = data,
               file = filename,
               ...)
    
    Sys.chmod(paths = filename,mode = "0777",use_umask = F)
  }
  
  # print("保存完毕")
}
