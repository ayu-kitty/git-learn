#!/opt/conda/bin/Rscript

#' 保存xlsx文件，对于check结果的标注
#'
#' @param data dataframe数据
#' @param filename 表格名
#' @param sheet sheet名
#'
#' @export
savexlsx8 <- function(data,
                      filename,
                      sheet = gsub(pattern = "\\.xlsx$",replacement = "",x = basename(filename)),
                      sheetmoudle = addsheet8,
                      overwrite = TRUE,
                      ...) {
  
  savexlsx(data = data,
           filename = filename,
           sheet = sheet,
           sheetmoudle = sheetmoudle,
           overwrite = overwrite,
           ...)
  
}

#' 在workbook中添加sheet，减少频繁保存表格问题，对于check结果的标注
#'
#' @param data dataframe数据
#' @param wb openxlsx的workbook类
#' @param sheet sheet名
#'
#' @export
addsheet8 <- function(data, 
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
  
  mode1 <- openxlsx::createStyle(fontName = "Arial", textDecoration = "bold")
  mode2 <- openxlsx::createStyle(fontName = "Arial")
  openxlsx::addStyle(wb, sheet = sheet, style = mode1, rows = 1, cols = 1:(dim(data)[2] + 1), gridExpand = T)
  openxlsx::addStyle(wb, sheet = sheet, style = mode1, rows = 2:(dim(data)[1] + 1), cols = 1, gridExpand = T)
  openxlsx::addStyle(wb, sheet = sheet, style = mode2, rows = 2:(dim(data)[1] + 1), cols = 2:(dim(data)[2] + 1), gridExpand = T)
  
  mode3 <- openxlsx::createStyle(fontName = "Arial", fgFill = "red")
  for (i in 1:dim(data)[1]) {
    for (j in 1:dim(data)[2]) {
      if (isFALSE(data[i, j])) {
        openxlsx::addStyle(wb, sheet = sheet, style = mode3, rows = i + 1, cols = j, gridExpand = T)
      }
    }
  }
  
  mode4 <- openxlsx::createStyle(fontName = "Arial", fgFill = "green")
  for (i in 1:dim(data)[1]) {
    for (j in 1:dim(data)[2]) {
      if (isTRUE(data[i, j])) {
        openxlsx::addStyle(wb, sheet = sheet, style = mode4, rows = i + 1, cols = j, gridExpand = T)
      }
    }
  }
  
  return(wb)
}
