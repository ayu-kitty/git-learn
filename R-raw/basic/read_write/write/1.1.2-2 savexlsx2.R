#!/opt/conda/bin/Rscript

#' 保存xlsx文件，无格式
#'
#' @param data dataframe数据
#' @param filename 表格名
#' @param sheet sheet名
#'
#' @export
savexlsx2 <- function(data,
                      filename,
                      sheet = gsub(pattern = "\\.xlsx$",replacement = "",x = basename(filename)),
                      sheetmoudle = addsheet2,
                      overwrite = TRUE,
                      ...) {
  
  savexlsx(data = data,
           filename = filename,
           sheet = sheet,
           sheetmoudle = sheetmoudle,
           overwrite = overwrite,
           ...)
  
}

#' 在workbook中添加sheet，减少频繁保存表格问题，无格式
#'
#' @param data dataframe数据
#' @param wb openxlsx的workbook类
#' @param sheet sheet名
#'
#' @export
addsheet2 <- function(data, 
                      wb = openxlsx::createWorkbook(),
                      sheet,
                      ...){
  
  if(is.null(data)) {return(wb)}
  if(dim(data)[1] == 0) {return(wb)}
  
  sheet <- fixsheetname(sheet = sheet)
  wb <- addwb(wb = wb,
              sheet = sheet,
              data = data,
              ...)
  
  return(wb)
}
