#!/opt/conda/bin/Rscript

#' @export
MetNormalizer <- function(input = "表达矩阵.xlsx",
                          saminfo = "批次信息.xlsx",
                          minqc = 0,
                          minsample = 0,
                          savepath = "批次矫正"){
  
  suppressWarnings(library(MetNormalizer))
  suppressWarnings(library(dplyr))
  suppressWarnings(library(tidyverse))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(reshape))
  if (!file.exists(savepath)) {
    dir.create(savepath)
  }
  
  afterpath <- "./校正后/"
  beforepath <- "./校正前/"
  
  rawdata <- readdata(input)
  
  spinfo <- readdata(saminfo)
  spinfo[,'class'] <- as.character(spinfo[,'class'])
  spinfo <- spinfo[,c('sample','order','class')]
  class <- select(spinfo,c('sample','class'))#%>%column_to_rownames("sample") #for PCA
  colnames(spinfo) <- c("sample.name","injection.order","class")
  spinfo$class <- ifelse(spinfo$class != "QC", "Subject", spinfo$class)
  savetxt(data = spinfo,filename = paste0(savepath,"/sample.info.csv"),sep = ",")
  
  sub_raw1 <- rawdata[c("ID","m/z","Retention time (min)")]
  colnames(sub_raw1) <- c("name","mz","rt")
  sub_raw2 <- rawdata[,spinfo$sample.name]
  sub_raw <- cbind(sub_raw1,sub_raw2)
  savetxt(data = sub_raw,paste0(savepath,"/data.csv"),sep = ",")
  savetxt(data = rawdata,filename = paste0(savepath,"/校正前/data.txt"))
  
  #---------
  setwd(savepath)
  metNor(ms1.data.name = "data.csv",
         sample.info.name = "sample.info.csv",
         minfrac.qc = minqc,
         minfrac.sample = minsample,
         optimization = TRUE,
         threads = 3)
  
  data_svr_normalization <- readdata(filename = "./svr_normalization_result/data_svr_normalization.csv",row.names = 1)
  # data_svr_normalization <- data_svr_normalization[,-1]
  data_svr_normalization <- rename(data_svr_normalization,c("QC.nor.rsd"="qc_rsd","sample.nor.rsd"="sample_rsd"))
  data_svr_normalization <- select(data_svr_normalization,-c("mz","rt"))
  datainfo <- rawdata[,!(colnames(rawdata) %in% spinfo$sample.name)]
  resultdata <- merge(datainfo,data_svr_normalization,by.x = "ID",by.y = "name")
  savetxt(data = resultdata,paste0("校正后/data_svr_normalization.csv"),sep = ",")

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
  svr_data <- select(data_svr_normalization,-c("name","qc_rsd","sample_rsd"))
  
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
  
  frequencydraw(file = data_svr_normalization,savepath = afterpath,imagetype = "png")
  
  result <- select(resultdata,-c("qc_rsd","sample_rsd"))
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
  parser$add_argument("-mq","--minqc",default = "./", help = "最小QC样本数")
  parser$add_argument("-ms","--minsample",default = "./", help = "最小样本数")
  
  args <- parser$parse_args()
  result <- do.call(what = MetNormalizer,args = args)
} 
