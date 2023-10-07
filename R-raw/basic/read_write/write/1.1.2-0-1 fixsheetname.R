#!/opt/conda/bin/Rscript

#' 修正31个字符以上的sheet名
#' 
#' @param sheet sheet名称
#'
#' @export
fixsheetname <- function(sheet){
  if (nchar(sheet) > 31) {
    maxn <- floor(nchar(sheet) / 3)
    n <- ceiling((nchar(sheet) - 31) / 3)
    sheet1 <- paste0(substr(sheet, 1, maxn - n), substr(sheet, maxn + 1, 2 * maxn - n), substr(sheet, 2 * maxn + 1, 3 * maxn - n))
  } else {
    sheet1 <- sheet
  }
  
  return(sheet1)
}
