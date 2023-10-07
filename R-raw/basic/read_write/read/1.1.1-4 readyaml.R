#!/opt/conda/bin/Rscript

#' 读取ymal格式的文件 
#' 
#' @param filename 文件名 
#' @param ... 见`readRDS`
#' 
#' @export 
readyaml <- function(filename,
                     stringsAsFactors = F,
                     row.names = NULL,
                     ...){
  tryCatch({
    if(file.exists(filename)){
     suppressMessages(library("yaml"))
      data <- read_yaml(file = filename,
                        ...)
      return(data)
    }else{ 
      stop(paste0(filename,"文件不存在"))
    }
  })
}

