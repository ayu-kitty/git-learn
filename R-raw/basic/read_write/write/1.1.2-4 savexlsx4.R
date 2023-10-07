#!/opt/conda/bin/Rscript

#' 保存xlsx文件，首行加粗，字体Arial,部分字体变红（对于通路富集结果）
#'
#' @param data dataframe数据
#' @param filename 表格名
#' @param sheet sheet名
#'
#' @export
savexlsx4 <- function(data,
                      filename,
                      sheet = gsub(pattern = "\\.xlsx$",replacement = "",x = basename(filename)),
                      sheetmoudle = addsheet4,
                      overwrite = TRUE,
                      ...) {
  
  savexlsx(data = data,
           filename = filename,
           sheet = sheet,
           sheetmoudle = sheetmoudle,
           overwrite = overwrite,
           ...)
  
}

#' 在workbook中添加sheet，减少频繁保存表格问题，首行加粗，字体Arial,部分字体变红（对于通路富集结果）
#'
#' @param data dataframe数据
#' @param wb openxlsx的workbook类
#' @param sheet sheet名
#'
#' @export
addsheet4 <- function(data, 
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
  mode3 <- openxlsx::createStyle(fontName = "Arial", fontColour = "red")
  openxlsx::addStyle(wb, sheet = sheet, style = mode1, rows = 1, cols = 1:dim(data)[2], gridExpand = T)
  openxlsx::addStyle(wb, sheet = sheet, style = mode2, rows = 2:(dim(data)[1] + 1), cols = 1:dim(data)[2], gridExpand = T)
  openxlsx::addStyle(wb, sheet = sheet, style = mode3, rows = which(is.na(data[, "URL"])) + 1, cols = 1:dim(data)[2], gridExpand = T)
  openxlsx::saveWorkbook(wb, file = name, overwrite = TRUE)
  
  return(wb)
}
