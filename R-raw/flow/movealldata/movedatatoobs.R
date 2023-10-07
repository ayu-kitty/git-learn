#!/opt/conda/bin/Rscript

#' movedatatoobs
#'
#' 上传obs
#'
#' @param data 项目信息
#' @param obsyear 存储位置
#' @param path 上传路径
#' @param ...
#'
#' @export
movedatatoobs<- function(data,
                         obsyear  = as.POSIXlt(Sys.Date())$year+1900,
                         path = "./raw/质谱数据/发送数据/") {
   system(paste0("source /etc/profile; tool_UpDataToObs -s 代谢实验部/",obsyear,
                " -m -p ",path," -f '",gsub("结题报告","原始数据",data$info$basic$项目报告),"'"),ignore.stdout = F,ignore.stderr = F,wait = F)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-d","--data",default = "../data", help = "data.RData")
  parser$add_argument("-o","--obsyear",default = "as.POSIXlt(Sys.Date())$year+1900", help = "存储位置")
  parser$add_argument("-p","--path",default = "./raw/质谱数据/发送数据/", help = "上传路径") 
  
  args <- parser$parse_args()
  
  result <- do.call(what = movedatatoobs,args = args) 
  
}
