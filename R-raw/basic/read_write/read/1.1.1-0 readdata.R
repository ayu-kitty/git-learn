#!/opt/conda/bin/Rscript

#' @export
readdata <- function(filename,
                     ...){
  
  if(!is.character(filename)){
    # warning("filename不是向量，直接返回原始值",immediate. = T)
    if(is.data.frame(filename)){
      return(filename)
    }else if(is.list(filename)){
      if(all(unlist(lapply(filename,is.character)))){
        if(all(unlist(lapply(filename,function(x){length(x) < 3})))){
          if(all(unlist(lapply(filename,function(x){grepl(pattern = "\\.",x = x[1])})))){
            data <- lapply(filename, readdata,...)
            return(data)
          }
          return(filename)
        }
        return(filename)
      }
      return(filename)
    }else{
      return(filename)
    }
  }
  
  data <- singlereaddata(filename,...)
  return(data)
  
}

#' 自动判断读取文件类型，读取相应文件
#' 
#' @param filename 文件名
#' @param sheet sheet序号或sheet名
#' @param filetype 限制文件类型
#' @param ... 见数据类型相关函数
#'
#' @export 
singlereaddata <- function(filename,
                           sheet = 1,
                           filetype = c("xlsx","xls","txt","csv","rds","yaml","imzML","mzML","RData"),
                           ...){
  tryCatch({
    
    if(!is.character(filename)){
      # warning("filename不是向量，直接返回原始值",immediate. = T)
      return(filename)
    }
    
    i <- 1
    
    if(!file.exists(filename[1])){
      stop(paste0(filename[1],"文件不存在"))
    }
    
    while(TRUE){
      if(grepl(pattern = paste0("\\.",filetype[i],"$"),x = filename[1])){
        break 
      }
      if(i > length(filetype)){
        stop(paste0(filename[1],"文件类型现不支持"))
      }
      i <- i + 1
    }
    
    if(grepl(pattern = "\\.txt$",x = filename[1])){
      
      txttry <- try({data <- readtxt(filename = filename[1],
                                     sep = "\t",
                                     ...)},silent = F)
      
      if("try-error" %in% class(txttry)){
        linedata <- readLines(filename[1],encoding = "utf-8")
        data <- data.frame(linedata[-1])
        row.names(data) <- linedata[-1]
      }
      
    }else if(grepl(pattern = "\\.csv$",x = filename[1])){
      
      txttry <- try({data <- readtxt(filename = filename[1],
                                     sep = ",",
                                     ...)},silent = F)
      
      if("try-error" %in% class(txttry)){
        linedata <- readLines(filename[1],encoding = "utf-8")
        data <- data.frame(linedata[-1])
        row.names(data) <- linedata[-1]
      }
      
    }else if(grepl(pattern = "\\.xls$",x = filename[1])|grepl(pattern = "\\.xlsx$",x = filename[1])){
      
      if(length(filename) > 1){
        data <- readxlsx(filename = filename[1],
                         sheet = filename[2],
                         ...)
      }else{
        data <- readxlsx(filename = filename[1],
                         sheet = sheet,
                         ...)
      }
      
    }else if(grepl(pattern = "\\.rds$",x = filename[1]) | grepl(pattern = "\\.RData$",x = filename[1])){
      
      data <- readrds(filename = filename[1],
                      ...)
      
    }else if(grepl(pattern = "\\.yaml$",x = filename[1])){
      
      data <- readyaml(filename = filename[1],
                       ...)
      
    }else if(grepl(pattern = "\\.imzML$",x = filename[1])){
      
      data <- readimzml(filename = filename,
                        ...)
      
    }else if(grepl(pattern = "\\.mzML$",x = filename[1])){
      
      data <- readmzml(filename = filename,
                       ...)
      
    }else{
      stop(paste0(filename,"文件类型现不支持"))
    }
    
    return(data)
    
  })
}
