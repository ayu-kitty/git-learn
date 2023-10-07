#!/opt/conda/bin/Rscript

#' 保存rds文件
#'
#' @param data 数据
#' @param filename 文件名
#'
#' @export
saverds <- function(data,
                    filename,
                    overwrite = T,
                    ...) {
  
  if(is.null(data)) {return()}
  
  print(paste0("在", filename, "中保存"))
  
  dirfilename <- createdir(filename = dirname(filename))

  if(!file.exists(filename) | overwrite){
    saveRDS(object = data,
            file = filename)
    
    Sys.chmod(paths = filename,mode = "0777",use_umask = F)
  }
  
  # print("保存完毕")
}
