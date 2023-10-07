#!/opt/conda/bin/Rscript

#' 读取txt格式的文件 
#' 
#' @param filename 文件名 
#' @param sep 列分割方式
#' @param fileEncoding 文件编码模式，默认为'utf-8'
#' @param stringsAsFactors 是否创建因子，默认为FALSE
#' @param check.names 是否检查列名，默认为FALSE
#' @param ... 见`read.delim`
#'
#' @return data.frame 
#' @export 
readtxt <- function(filename,
                    sep = "\t",
                    fileEncoding = "UTF-8",
                    stringsAsFactors = FALSE,
                    check.names = FALSE,
                    row.names = NULL,
                    nrows = -1,
                    colClasses = NA,
                    ...){
  tryCatch({
    if(file.exists(filename)){
      if(!is.null(row.names)){
        
        data2 <- read.delim(file = filename,
                            sep = sep,
                            fileEncoding = fileEncoding,
                            stringsAsFactors = stringsAsFactors,
                            check.names = check.names,
                            encoding = "UTF-8", 
                            row.names = NULL,
                            nrows = 1,
                            colClasses = colClasses,
                            ...)
        
        colClasses <- rep(x = NA,dim(data2)[2])
        colClasses[row.names] <- "character"
        
        data <- read.delim(file = filename,
                           sep = sep,
                           fileEncoding = fileEncoding,
                           stringsAsFactors = stringsAsFactors,
                           check.names = check.names,
                           encoding = "UTF-8", 
                           row.names = row.names,
                           colClasses = colClasses,
                           nrows = nrows,
                           ...)
        
      }else{
        data <- read.delim(file = filename,
                           sep = sep,
                           fileEncoding = fileEncoding,
                           stringsAsFactors = stringsAsFactors,
                           check.names = check.names,
                           encoding = "UTF-8", 
                           row.names = row.names,
                           nrows = nrows,
                           colClasses = colClasses,
                           ...)
      }
      return(data)
    }else{ 
      stop(paste0(filename,"文件不存在"))
    }
  })
}
