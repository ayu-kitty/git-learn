#!/opt/conda/bin/Rscript

#' 保存xlsx文件，首行加粗，字体Arial
#'
#' @param data dataframe数据
#' @param filename 表格名
#' @param sheet sheet名
#'
#' @export
savexlsx <- function(data,
                     filename,
                     sheet = gsub(pattern = "\\.xlsx$",replacement = "",x = basename(filename)),
                     sheetmoudle = addsheet1,
                     overwrite = TRUE,
                     ...) {
  if(is.null(data)) {return()}
  
  print(paste0("在", filename, "中保存", sheet, "表格"))
  
  if (file.exists(filename)) {wb <- openxlsx::loadWorkbook(filename)
  } else {wb <- openxlsx::createWorkbook()}
  
  wb <- addsheet(data = data,
                 wb = wb,
                 sheet = sheet,
                 sheetmoudle = sheetmoudle,
                 ...)
  
  savewb(wb = wb,filename = filename,overwrite = overwrite)
  
  # print("保存完毕")
}

#' @export
addsheet <- function(data,
                     wb,
                     sheet = NULL,
                     sheetmoudle = addsheet1,
                     ...){
  if(is.data.frame(data)){
    
    if(is.null(sheet)){sheet <- "sheet1"}
    
    wb <- sheetmoudle(data = data,
                      wb = wb,
                      sheet = sheet,
                      ...)
    
  }else if(is.list(data)){
    
    for ( i in 1:length(data)) {
      
      if((is.null(names(data)))|("" %in% names(data)[i])){
        if(is.null(sheet)){
          sheet1 <- "sheet"
          sheet1 <- paste0(c(sheet1,i),collapse = "")
        }else{
          sheet1 <- sheet
          sheet1 <- paste0(c(sheet1,i),collapse = "-")
        }
      }else{
        sheet1 <- names(data)[i]
      }
      
      wb <- addsheet(data = data[[i]],
                     wb = wb,
                     sheet = sheet1,
                     sheetmoudle = sheetmoudle,
                     ...)
    }
  }
  
  return(wb)
}


#' @export
addwb <- function(wb,
                  sheet,
                  data,
                  ...){
  
  tryfetch <- try({
    openxlsx::addWorksheet(wb, sheet)
  },silent = T)
  
  if ("try-error" %in% class(tryfetch)) {
    openxlsx::removeWorksheet(wb, sheet)
    openxlsx::addWorksheet(wb, sheet)
  }
  
  openxlsx::writeData(wb = wb, 
                      sheet = sheet, 
                      x = data,
                      ...)
  
  return(wb)
}

#' 保存workbook对象
#' 
#' @param wb workbook对象
#' @param filename 保存xlsx文件名
#'
#' @export
savewb <- function(wb,filename,overwrite = TRUE){
  
  dirfilename <- createdir(filename = dirname(filename))
  
  if(!file.exists(filename) | overwrite){
    openxlsx::saveWorkbook(wb = wb, 
                           file = filename, 
                           overwrite = overwrite)
  }
  
  Sys.chmod(paths = filename,mode = "0777",use_umask = F)
}