#!/opt/conda/bin/Rscript

#' 保存xlsx文件，首行加粗，字体Arial,上调红色，下调蓝色
#'
#' @param data dataframe数据
#' @param name 表格名
#' @param sheet sheet名
#'
#' @export
savexlsx_fc <- function(data,
                        filename,
                        sheet = gsub(pattern = "\\.xlsx$",replacement = "",x = basename(filename)),
                        sheetmoudle = addsheet_fc,
                        overwrite = TRUE,
                        ...) {
  
  savexlsx(data = data,
           filename = filename,
           sheet = sheet,
           sheetmoudle = sheetmoudle,
           overwrite = overwrite,
           ...)
  
}

#' 在workbook中添加sheet，减少频繁保存表格问题
#'
#' @param data dataframe数据
#' @param wb openxlsx的workbook类
#' @param sheet sheet名
#'
#' @export
addsheet_fc <- function(data, 
                        wb,
                        sheet,
                        ...){
  
  if(is.null(data)) {return(wb)}
  if(dim(data)[1] == 0) {return(wb)}
  
  sheet1 <- fixsheetname(sheet = sheet)
  
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
  
  wb1 <- addwb(wb = wb,
               sheet = sheet1,
               data = data,
               ...)
  
  mode1 <- openxlsx::createStyle(fontName = "Arial", textDecoration = "bold",halign = "left")
  mode2 <- openxlsx::createStyle(fontName = "Arial",halign = "left")
  # mode3 <- openxlsx::createStyle(fontName = "Arial", fgFill = "powderblue")
  # mode4 <- openxlsx::createStyle(fontName = "Arial", fgFill = "brown1")
  mode3 <- openxlsx::createStyle(fontName = "Arial", fgFill = "#727fb5")
  mode4 <- openxlsx::createStyle(fontName = "Arial", fgFill = "#f47d8c")
  openxlsx::addStyle(wb1, sheet = sheet1, style = mode1, rows = 1, cols = 1:dim(data)[2], gridExpand = T)
  openxlsx::addStyle(wb1, sheet = sheet1, style = mode2, rows = 2:(dim(data)[1] + 1), cols = 1:dim(data)[2], gridExpand = T)
  
  if("log2(FC)" %in% colnames(data)){
    openxlsx::addStyle(wb1, sheet = sheet1, style = mode3, rows = which(data[, "log2(FC)"] < 0) + 1, cols = which(colnames(data) == "log2(FC)"), gridExpand = T)
    openxlsx::addStyle(wb1, sheet = sheet1, style = mode4, rows = which(data[, "log2(FC)"] > 0) + 1, cols = which(colnames(data) == "log2(FC)"), gridExpand = T)
  }
  if("FC" %in% colnames(data)){
    openxlsx::addStyle(wb1, sheet = sheet1, style = mode3, rows = which(data[, "FC"] < 1) + 1, cols = which(colnames(data) == "FC"), gridExpand = T)
    openxlsx::addStyle(wb1, sheet = sheet1, style = mode4, rows = which(data[, "FC"] > 1) + 1, cols = which(colnames(data) == "FC"), gridExpand = T)
  }
  if("log2FoldChange" %in% colnames(data)){
    openxlsx::addStyle(wb1, sheet = sheet1, style = mode3, rows = which(data[, "log2FoldChange"] < 0) + 1, cols = which(colnames(data) == "log2FoldChange"), gridExpand = T)
    openxlsx::addStyle(wb1, sheet = sheet1, style = mode4, rows = which(data[, "log2FoldChange"] > 0) + 1, cols = which(colnames(data) == "log2FoldChange"), gridExpand = T)
  }
  if("FoldChange" %in% colnames(data)){
    openxlsx::addStyle(wb1, sheet = sheet1, style = mode3, rows = which(data[, "FoldChange"] < 1) + 1, cols = which(colnames(data) == "FoldChange"), gridExpand = T)
    openxlsx::addStyle(wb1, sheet = sheet1, style = mode4, rows = which(data[, "FoldChange"] > 1) + 1, cols = which(colnames(data) == "FoldChange"), gridExpand = T)
  }
  if("Regulation" %in% colnames(data)){
    openxlsx::addStyle(wb1, sheet = sheet1, style = mode3, rows = which(data[, "Regulation"] == "Down") + 1, cols = which(colnames(data) == "Regulation"), gridExpand = T)
    openxlsx::addStyle(wb1, sheet = sheet1, style = mode4, rows = which(data[, "Regulation"] == "Up") + 1, cols = which(colnames(data) == "Regulation"), gridExpand = T)
  }
  
  if("cidlink" %in% colnames(data)){
    linkStyle <- openxlsx::createStyle(fontColour = "#0000FF", textDecoration = "underline")
    openxlsx::addStyle(wb, sheet= sheet1, style = linkStyle, rows = 2:(dim(data)[1] + 1), cols = which(colnames(data)=="cidlink"))
    openxlsx::setColWidths(wb, sheet= sheet1, widths = 45,cols = which(colnames(data)=="cidlink"))
  }
  
  return(wb1)
}

