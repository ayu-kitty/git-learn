#!/opt/conda/bin/Rscript

#' 空带Qualitative中3个class计数
#' 
#' @param filepath 文件读取路径
#' @param sheet 对于Qualitative.xlsx中需要处理的sheet名
#' 
#' @export
classnumsp_dataextraction <- function(filepath = "./",
                                      sheet = "pos-all"){
  data <- read.xlsx(paste0(filepath,"Qualitative.xlsx"),sheet = sheet)
  superclass <- unique(data$Super.Class)
  subclass <- unique(data$Sub.Class)
  class <- unique(data$Class)
  
  data0 <- list()
  data0[["superclass"]] <- data.frame(superclass = superclass,
                                      superclass_num = unlist(lapply(superclass, function(x){length(grep(x,data$Super.Class))})),
                                      superclass_percent = paste0(unlist(lapply(superclass, function(x){round(100*length(grep(x,data$Super.Class))/nrow(data),digits = 2)})),"%"))
  data0[["subclass"]] <- data.frame(subclass = subclass,
                                    subclass_num = unlist(lapply(subclass, function(x){length(grep(x,data$Sub.Class))})),
                                    subclass_percent = paste0(unlist(lapply(subclass, function(x){round(100*length(grep(x,data$Sub.Class))/nrow(data),digits = 2)})),"%"))
  data0[["class"]] <- data.frame(class = class,
                                 class_num = unlist(lapply(class, function(x){length(grep(x,data$Class))})),
                                 class_percent = paste0(unlist(lapply(class, function(x){round(100*length(grep(x,data$Class))/nrow(data),digits = 2)})),"%"))
  lmbio::savexlsx(data0,filename = paste0(filepath,"代谢物类别信息.xlsx"),sheet = sheet)
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-wp","--filepath",default = "./", help = "文件读取路径")
  parser$add_argument("-m","--sheet",default = "pos-all", help = "对于Qualitative.xlsx中需要处理的sheet名")
 
  args <- parser$parse_args()
  
  writeinfo()
  
  mulargs <- do.call(what = classnumsp_dataextraction,args = args)
  
  writeinfo(endtime = T)
  
}
