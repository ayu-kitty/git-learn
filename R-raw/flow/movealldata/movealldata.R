#!/opt/conda/bin/Rscript

#' movealldata
#'
#' 上传原始数据
#'
#' @param data 项目信息
#' @param ...
#'
#' @export

movealldata <- function(data,
                           ...) {
  UseMethod("movealldata")
}

movealldata.default<-function(data,
                      ...){
  
  load("raw.RData")
  if(data[["info"]][["basic"]][["S3"]] %in% c("UntargetBothPmdb","UntargetBothEmdb","UntargetBoth")){
  #将原始数据移动到发送文件夹中 
  movedata(lcdata)
  #上机对照表整理并放在发送数据文件夹中
  if(length(list.files("./raw",pattern = "上机对照表")) == 1){
    modeldata <- readxlsx(filename = "./raw/上机对照表.xlsx",sheet = 1)
    mzmlsamplename <- lcdata$info$sample$samplename[!(lcdata$info$sample$deal %in% "删除")]
    modeldata2 <- modeldata[mzmlsamplename %in% modeldata$样本分析名,]
    lmbio::savexlsx1(modeldata2,filename = "./raw/质谱数据/发送数据/上机对照表.xlsx",sheet = "样本对照表")
    }
    #将原始数据移动到发送文件夹中 
  movedata(gcdata)
  #上机对照表整理并放在发送数据文件夹中
  if(length(list.files("./raw",pattern = "上机对照表")) == 1){
    modeldata <- readxlsx(filename = "./raw/上机对照表.xlsx",sheet = 1)
    mzmlsamplename <- gcdata$info$sample$samplename[!(gcdata$info$sample$deal %in% "删除")]
    modeldata2 <- modeldata[mzmlsamplename %in% modeldata$样本分析名,]
    lmbio::savexlsx1(modeldata2,filename = "./raw/质谱数据/发送数据/上机对照表.xlsx",sheet = "样本对照表")
  }
}else{
  #将原始数据移动到发送文件夹中 
  movedata(data)
  #上机对照表整理并放在发送数据文件夹中
  if(length(list.files("./raw",pattern = "上机对照表")) == 1){
    modeldata <- readxlsx(filename = "./raw/上机对照表.xlsx",sheet = 1)
    mzmlsamplename <- data$info$sample$samplename[!(data$info$sample$deal %in% "删除")]
    modeldata2 <- modeldata[mzmlsamplename %in% modeldata$样本分析名,]
    lmbio::savexlsx1(modeldata2,filename = "./raw/质谱数据/发送数据/上机对照表.xlsx",sheet = "样本对照表")
  }
  }
  #原始数据md5验证及上传obs 
  movedatatoobs(data) 

}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-d","--data",default = "../data", help = "data.RData")
  args <- parser$parse_args()
  
  result <- do.call(what = movealldata,args = args) 
  
}

