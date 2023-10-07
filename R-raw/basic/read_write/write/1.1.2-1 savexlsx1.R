#!/opt/conda/bin/Rscript

#' 保存xlsx文件，首行加粗，字体Arial
#'
#' @param data dataframe数据
#' @param filename 表格名
#' @param sheet sheet名
#'
#' @export
savexlsx1 <- function(data,
                      filename,
                      sheet = gsub(pattern = "\\.xlsx$",replacement = "",x = basename(filename)),
                      sheetmoudle = addsheet1,
                      overwrite = TRUE,
                      ...) {
  
  savexlsx(data = data,
           filename = filename,
           sheet = sheet,
           sheetmoudle = sheetmoudle,
           overwrite = overwrite,
           ...)
  
}

#' 在workbook中添加sheet，减少频繁保存表格问题，首行加粗，字体Arial
#'
#' @param data dataframe数据
#' @param wb openxlsx的workbook类
#' @param sheet sheet名
#'
#' @export
addsheet1 <- function(data, 
                      wb = openxlsx::createWorkbook(),
                      sheet,
                      cid = T,
                      ...){
  
  if(is.null(data)) {return(wb)}
  if(dim(data)[1] == 0) {return(wb)}
  
  sheet <- fixsheetname(sheet = sheet)
  if(cid){
    if("cid" %in% colnames(data)){
      if(any(grepl(pattern = ";",data$cid))){
        
      }else{
        if(!("cidlink" %in% colnames(data))){
          data1 <- data[,1:which(colnames(data)=="cid"),drop= F]
          data1[,"cidlink"] <- NA
          if(which(colnames(data)=="cid") != dim(data)[2]){
            data2 <- data[,(which(colnames(data)=="cid")+1):dim(data)[2],drop= F]
            data <- cbind(data1,data2)
          }else{
            data <- data1
          }
        }
        
        for ( i in 1:dim(data)[1]) {
          if(!is.na(data$cid[i])){
            data$cidlink[i] <- paste0("HYPERLINK(\"",
                                      "http://funmeta.oebiotech.com/network/",data$cid[i],
                                      "\", \"",
                                      "http://funmeta.oebiotech.com/network/",data$cid[i],
                                      "\")")
          }
        }
        
        class(data$cidlink) <- "formula"
      }
    }
  }
  
  # base::print(data)
  
  wb <- addwb(wb = wb,
              sheet = sheet,
              data = data,
              ...)
  
  mode1 <- openxlsx::createStyle(fontName = "Arial", textDecoration = "bold",halign = "left")
  mode2 <- openxlsx::createStyle(fontName = "Arial",halign = "left")
  openxlsx::addStyle(wb, sheet = sheet, style = mode1, rows = 1, cols = 1:dim(data)[2], gridExpand = T)
  openxlsx::addStyle(wb, sheet = sheet, style = mode2, rows = 2:(dim(data)[1] + 1), cols = 1:dim(data)[2], gridExpand = T)
  
  if("cidlink" %in% colnames(data)){
    linkStyle <- openxlsx::createStyle(fontColour = "#0000FF", textDecoration = "underline")
    openxlsx::addStyle(wb, sheet= sheet, style = linkStyle, rows = 2:(dim(data)[1] + 1), cols = which(colnames(data)=="cidlink"))
    openxlsx::setColWidths(wb, sheet= sheet, widths = 45,cols = which(colnames(data)=="cidlink"))
  }
  
  return(wb)
}
