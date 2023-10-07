#!/opt/conda/bin/Rscript

#' 读取xlsx、xls格式的文件
#' 
#' @param filename 文件名 
#' @param sheet sheet序号或sheet名
#' @param stringsAsFactors 是否创建因子，默认为FALSE
#' @param row.names 使用第几列当行名，默认为NULL，无行名
#' @param ... 见`readxl::read_excel`
#'
#' @return data.frame 
#' @export 
readxlsx <- function(filename,
                     sheet = 1,
                     stringsAsFactors = F,
                     row.names = NULL,
                     ...){
  suppressMessages(library("readxl"))
  
  if(is.list(sheet)){
    stop("sheet不能输入list")
  }

  if(length(sheet) > 1){
    data <- lapply(sheet, function(sheet,...){readxlsx(sheet = sheet,...)},
                   filename = filename,
                   stringsAsFactors = stringsAsFactors,
                   row.names = row.names,
                   ...)
    names(data) <- getsheetname(filename = filename,sheet = sheet)
    return(data)
  }
  
  tryCatch({
    if(file.exists(filename)){
      
      sheet <- getsheetname(filename = filename,sheet = sheet)
      
      if(sheet == ""){
        data <- readtxt(filename = filename,
                        stringsAsFactors = stringsAsFactors,
                        row.names = row.names,
                        ...)
      }else{
        data <- read_excel(path = filename,
                           sheet = sheet,
                           progress = F,
                           guess_max = 20000,
                           ...)
        
        data <- as.data.frame(data,
                              stringsAsFactors = stringsAsFactors)
        if(!is.null(row.names)){
          if(is.logical(row.names)){
            if(row.names){
              row.names <- 1
            }else{
              return(data)
            }
          }
          if(!any(is.na(data[,row.names]))){
            if(!any(duplicated(data[,row.names]))){
              row.names(data) <- data[,row.names]
              data <- data[,-row.names,drop = F]
            }else{
              stop(paste0("在",filename,"文件中",sheet,"表第",row.names,"列有重复不能作为行名，请检查"))
            }
          }else{
            stop(paste0("在",filename,"文件中",sheet,"表第",row.names,"列有空值，请检查"))
          }
        }
      }
      
      return(data)
    }else{ 
      stop(paste0(filename,"文件不存在"))
    }
  })
}


#' 获取sheet名
#'
#' @param filename 文件名 
#' @param sheet sheet名称
#'
#' @return sheet名向量
#' 
#' @export
getsheetname <- function(filename,
                         sheet = NULL) {
  suppressMessages(library("readxl"))
  
  tryCatch({
    
    if(!is.vector(filename)){
      # warning("filename不是向量，直接返回空",immediate. = T)
      return("")
    }
    
    if(length(filename) > 1){
      filename <- filename[1]
    }
    
    if(is.list(sheet)){
      sheetname <- lapply(sheet, function(sheet,...){getsheetname(sheet =sheet,...)},
                          filename = filename)
      return(sheetname)
    }
    
    if(grepl(pattern = "\\.xls$",x = filename)|grepl(pattern = "\\.xlsx$",x = filename)){
      if(file.exists(filename)){
        tryinfo <- try({
          sheetname <- excel_sheets(filename)
        },silent = T)
        
        if("try-error" %in% class(tryinfo)){
          sheetname <- ""
          return(sheetname)
        }
        
        if(is.null(sheet)){
          
        }else if(is.numeric(sheet)){
          if(all(sheet %in% 1:length(sheetname))){
            sheetname <- sheetname[sheet]
          }else{
            stop(paste0("在",filename,"文件中sheet选择超出范围"))
          }
        }else{
          if(all(sheet %in% sheetname)){
            sheetname <- sheet
          }else{
            sheet <- sheet[!(sheet %in% sheetname)]
            
            stop(paste0("在",filename,"文件中无以下sheet:",paste(sheet,collapse = ";")))
          }
        }
        
      }else{
        stop(paste0(filename,"文件不存在"))
      }
    }else{
      sheetname <- ""
    }
    
  })
  
  return(sheetname)
}

