#!/opt/conda/bin/Rscript

#' @export
Stattarget <- function(input = "表达矩阵.xlsx",
                       saminfo = "批次信息.xlsx",
                       savepath = "批次矫正",
                       Frule = 0.8, 
                       MLmethod = "QCRFSC", 
                       QCspan = 0,
                       imputeM = "KNN"){
  
  suppressWarnings(library(statTarget))
  # suppressWarnings(library(tidyverse))
  suppressWarnings(library(dplyr))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(reshape))
  
  if (!file.exists(savepath)) {
    dir.create(savepath)
  }
  afterpath <- "./校正后/"
  beforepath <- "./校正前/"
  
  rawdata <- readdata(filename=input)
  
  sip <- readdata(saminfo)
  class <- sip
  sip$class <- gsub("QC","NA",sip$class)
  sip <<- sip
  savetxt(data = sip,paste0(savepath,"/sample.info.csv"),sep = ",")
  
  sub_raw1 <- rawdata[,c("ID"),drop = F]
  colnames(sub_raw1) <- c("name")
  sub_raw2 <- rawdata[,sip$sample]
  sub_raw <- cbind(sub_raw1,sub_raw2)
  # sub_raw  <- sub_raw[!apply(sub_raw[,-1,drop = F],MARGIN = 1,FUN = function(x){all(x == 0)}),]
  datainfo <- rawdata[,!(colnames(rawdata) %in% sip$sample)]
  savetxt(data = sub_raw,paste0(savepath,"/data.csv"),sep = ",")
  savetxt(data = rawdata,filename = paste0(savepath,"/校正前/data.txt"))
  
  #---------批次矫正
  setwd(savepath)
  shiftCor(samPeno = "sample.info.csv",
           samFile = "data.csv", 
           Frule = Frule,
           MLmethod = MLmethod, 
           QCspan = QCspan,
           imputeM = imputeM) 
  #---------
  
  # afterpath <- "statTarget/shiftCor/After_shiftCor/"
  # beforepath <- "statTarget/shiftCor/Before_shiftCor/"
  correctresult <- read.csv(file = "statTarget/shiftCor/After_shiftCor/shift_all_cor.csv",row.names = "sample",header=T,check.names=F)
  svr_data <- t(correctresult[,-1])
  svr_data <- as.data.frame(svr_data)
  result <- merge(datainfo,svr_data,by.x = "ID",by.y = 0)
  
  #将分组文件转化为list,为PCA class输入文件
  my_class <- split(class, class$class)
  my_class <- lapply(my_class, function(x) x$sample)
  
  #PCA原始数据
  raw_pca <- mulstatistics_file(datafile = sub_raw2,
                                classfile = my_class,
                                group = unique(class$class),
                                mode = "PCA")
  
  map_mulstatistics_scoremap(filename = raw_pca,
                             mapname = "raw_pca",
                             savepath = beforepath,
                             imagetype = "png")
  
  map_common_corrplot2(filename = t(sub_raw2),
                       mapname = "raw_corrplot",
                       savepath = beforepath,
                       imagetype = "png")
  
  sub_raw2[,"qc_rsd"] <- apply(sub_raw2[,class$sample[class$class == "QC"],drop = F],1,function(x){sd(x)/mean(x)*100})
  frequencydraw(file = sub_raw2,savepath = beforepath,imagetype = "png")
  
  #PCA矫正后
  svr_pca <- mulstatistics_file(datafile = svr_data,
                                classfile = my_class,
                                mode = "PCA", 
                                group = unique(class$class))
  
  map_mulstatistics_scoremap(filename = svr_pca,
                             mapname = "svr_pca",
                             savepath = afterpath,
                             imagetype = "png")
  
  map_common_corrplot2(filename = t(svr_data),
                       mapname = "srv_corrplot",
                       savepath = afterpath,
                       imagetype = "png")
  
  svr_data[,"qc_rsd"] <- apply(svr_data[,class$sample[class$class == "QC"],drop = F],1,function(x){sd(x)/mean(x)*100})
  frequencydraw(file = svr_data,savepath = afterpath,imagetype = "png")
  
  savetxt(data = result,filename = paste0(afterpath,"/data.txt"))
  system("chmod -R 777 ./* ")
  return(result)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-i","--input",default = "表达矩阵.xlsx", help = "表达矩阵")
  parser$add_argument("-si","--saminfo",default = "批次信息.xlsx", help = "样本批次信息数据，包括样本名称、批次、分类和上机顺序")
  parser$add_argument("-sp","--savepath",default = "批次矫正", help = "保存数据路径")
  parser$add_argument("-fr","--Frule",default = 0.8, help = "80%过滤原则")
  parser$add_argument("-m","--MLmethod",default = "QCRFSC", help = "基于QC的信号校正的机器学习方法,QCRFSC:随机森林,QCRLSC:LOESS信号")
  parser$add_argument("-qs","--QCspan",default = 0, help = "QCRLSC的平滑参数，用于控制QCRLSS方法中的偏差-方差平衡,为0表示广义交叉验证,防止过拟合")
  parser$add_argument("-im","--imputeM",default = "KNN", help ="缺失值填充方法")
  
  args <- parser$parse_args()
  result <- do.call(what = Stattarget,args = args) 
} 
