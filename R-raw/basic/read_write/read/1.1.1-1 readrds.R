#!/opt/conda/bin/Rscript

#' 读取rds格式的文件 
#' 
#' @param filename 文件名 
#' @param ... 见`readRDS`
#' 
#' @export 
readrds <- function(filename,
                    stringsAsFactors = F,
                    row.names = NULL,
                    ...){
  tryCatch({
    if(file.exists(filename)){
      data <- readRDS(file = filename,
                      ...)
      return(data)
    }else{ 
      stop(paste0(filename,"文件不存在"))
    }
  })
}
