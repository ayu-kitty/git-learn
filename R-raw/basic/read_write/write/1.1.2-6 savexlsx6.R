#!/opt/conda/bin/Rscript

#' 保存xlsx文件，通路富集结果中代谢物标注颜色使用
#'
#' @param data dataframe数据
#' @param filename 表格名
#' @param sheet sheet名
#'
#' @export
savexlsx6 <- function(data,
                      filename,
                      sheet = gsub(pattern = "\\.xlsx$",replacement = "",x = basename(filename)),
                      sheetmoudle = addsheet6,
                      overwrite = TRUE,
                      ...) {
  
  savexlsx(data = data,
           filename = filename,
           sheet = sheet,
           sheetmoudle = sheetmoudle,
           overwrite = overwrite,
           ...)
  
}

#' 在workbook中添加sheet，减少频繁保存表格问题，通路富集结果中代谢物标注颜色使用
#'
#' @param data dataframe数据
#' @param wb openxlsx的workbook类
#' @param sheet sheet名
#'
#' @export
addsheet6 <- function(data, 
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
  
  openxlsx::writeData(wb, sheet = sheet, data)
  openxlsx::writeData(wb, sheet = "说明", data.frame("说明" = c("绿色背景为富集到通路的差异代谢物", "蓝色背景为有KEGG ID，但没有富集到通路的差异代谢物", "无背景色为没有KEGG ID的差异代谢物")))
  mode1 <- openxlsx::createStyle(fontName = "Arial", textDecoration = "bold")
  mode2 <- openxlsx::createStyle(fontName = "Arial")
  mode3 <- openxlsx::createStyle(fontName = "Arial", fgFill = "palegreen")
  mode4 <- openxlsx::createStyle(fontName = "Arial", fgFill = "powderblue")
  openxlsx::addStyle(wb, sheet = sheet, style = mode1, rows = 1, cols = 1:dim(data)[2], gridExpand = T)
  openxlsx::addStyle(wb, sheet = sheet, style = mode2, rows = 2:(dim(data)[1] + 1), cols = 1:dim(data)[2], gridExpand = T)
  openxlsx::addStyle(wb, sheet = sheet, style = mode4, rows = c(2:(dim(data)[1] + 1))[!is.na(data[, "KEGG"])], cols = 1:dim(data)[2], gridExpand = T)
  openxlsx::addStyle(wb, sheet = sheet, style = mode3, rows = c(2:(dim(data)[1] + 1))[!is.na(data[, "ID Annotation"])], cols = 1:dim(data)[2], gridExpand = T)
  
  openxlsx::addStyle(wb, sheet = "说明", style = mode1, rows = 1, cols = 1, gridExpand = T)
  openxlsx::addStyle(wb, sheet = "说明", style = mode2, rows = 4, cols = 1, gridExpand = T)
  openxlsx::addStyle(wb, sheet = "说明", style = mode4, rows = 3, cols = 1, gridExpand = T)
  openxlsx::addStyle(wb, sheet = "说明", style = mode3, rows = 2, cols = 1, gridExpand = T)
  
  return(wb)
}
