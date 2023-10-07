#!/opt/conda/bin/Rscript

#' @export
map_autodraw <- R6::R6Class("map_autodraw",
                            public = list(
                              initialize = function(moudle,row.names = NULL){
                                self$moudle <- moudle
                                self$row.names <- row.names
                              },
                              draw = function(...){
                                autodraw(moudle = self$moudle,
                                         row.names = self$row.names,
                                         ...)
                              },
                              moudle = NULL,
                              row.names = NULL
                            )
)


#' 根据提供文件及函数进行批量绘图
#' 
#' @param filename 文件名 
#' @param sheet sheet名
#' @param moudle 绘图函数
#' @param errorstop 逻辑，是否报错停止
#' @param frontname 前置名称
#' @param backname 后置名称
#' @param row.names 使用第几列当行名，默认为NULL，无行名
#' @param ... 见moudle函数
#'
#' @export
autodraw <- function(filename,
                     ...) {
  
  if(is.character(filename)){
    args <- lapply(X = filename,
                   FUN = autosingledraw,
                   ...)
  }else{
    args <- autosingledraw(filename = filename,...)
  }
  
}


#' 单个文件进行可视化
#'
#' @param filename 文件名 
#' @param sheet sheet名
#' @param moudle 绘图函数
#' @param errorstop 逻辑，是否报错停止
#' @param frontname 前置名称
#' @param backname 后置名称
#' @param row.names 使用第几列当行名，默认为NULL，无行名
#' @param ... 见moudle函数
#'
#' @export
autosingledraw <- function(filename,
                           sheet = NULL,
                           moudle = NULL,
                           errorstop = T,
                           mapname = NULL,
                           row.names = NULL,
                           ...) {
  
  if (!is.function(moudle)) {
    stop("moudle模块未导入成功")
  }
  
  sheetname <- getsheetname(filename = filename,
                            sheet = sheet)
  
  for (i in 1:length(sheetname)) {
    plottry <- try({
      if(is.character(filename)){
        print(paste0(filename,"文件",sheetname[[i]],"运行开始"))
      }else{
        print("运行开始")
      }
      
      data <- readdata(filename = filename,
                       sheet = sheetname[[i]],
                       row.names = row.names)
      
      if(is.null(mapname)){
        if(is.character(filename)){
          mapname2 <- gsub(pattern = "\\..*$",replacement = "",x = basename(filename))
        }else{
          mapname2 <- NULL
        }
      }else{
        mapname2 <- mapname
      }
      
      if(length(sheetname) > 1 & is.character(sheetname)){
        backname <- paste(sheetname[[i]],collapse = "_")
      }else if("" %in% sheetname[[i]]){
        backname <- NULL
      }else{
        backname <- NULL
      }
      
      if(is.character(filename)){
        returndata <- moudle(data = data,
                             mapname = paste(c(mapname2, backname),collapse = "-"),
                             ...)
      }else{
        if(is.null(mapname)){
          returndata <- moudle(data = data,
                               ...)
        }else{
          returndata <- moudle(data = data,
                               mapname = mapname,
                               ...)
        }
      }
      
      # if(is.character(filename)){
      #   print(paste0(filename,"文件",sheetname[[i]],"运行完成"))
      # }else{
      #   print("运行完成")
      # }

    },silent = F)
    
    if ("try-error" %in% class(plottry)& errorstop) {
      if(is.character(filename)){
        stop(paste0(filename,"文件",sheetname[[i]],"运行失败"))
      }else{
        stop("运行失败")
      }
    } else if ("try-error" %in% class(plottry)) {
      if(is.character(filename)){
        warning(paste0(filename,"文件",sheetname[[i]],"运行失败"), immediate. = T)
      }else{
        warning("运行失败", immediate. = T)
      }
    }
  }
  
  return("运行完成")
}
