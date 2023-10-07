#!/opt/conda/bin/Rscript
#' @export
SERRF_run <- function(input = "表达矩阵.xlsx",
                      saminfo = "批次信息.xlsx",
                      savepath = "批次矫正",
                      ...){

  suppressWarnings(library(SERRFweb))
  suppressWarnings(library(dplyr))
  suppressWarnings(library(tidyverse))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(reshape))
  if (!file.exists(savepath)) {
    dir.create(savepath)
  }
  afterpath <- paste0(savepath,"/校正后/")
  beforepath <- paste0(savepath,"/校正前/")
  
  
  #将数据矩阵和样本信息合并
  rawdata <- readdata(filename = input)
  sub_raw1 <- rawdata[,0]
  sub_raw1[,"No"] <- 1:dim(sub_raw1)
  sub_raw1[,"label"] <- rawdata[,"ID"]
  
  #样本信息
  SP <- readdata(saminfo)
  sip <- SP
  class <- sip
  colnames(SP) <- c("label","batch","sampleType","time")
  class <- select(SP,c("label","sampleType")) #for PCA
  SP$sampleType <- ifelse(SP$sampleType != "QC", "sample", SP$sampleType)
  SP$sampleType <- gsub("QC","qc",SP$sampleType)
  SP <- rbind(colnames(SP),SP)
  SP <- SP[,c("batch","sampleType","time","label")]
  df_transposed <- t(SP)
  colnames(df_transposed) <- df_transposed[4,]
  df_transposed <- cbind("No" = rep("", nrow(df_transposed)), df_transposed)
  
  sub_raw2 <- rawdata[,sip$sample]
  sub_raw <- cbind(sub_raw1,sub_raw2)
  
  all_data <- rbind(df_transposed,sub_raw)
  all_data[4,1] <- "No" 
  
  savetxt(data = rawdata,filename = paste0(savepath,"/校正前/data.txt"))
  openxlsx::write.xlsx(x = all_data,file = paste0(savepath,"/all_data.xlsx"),rowNames = FALSE,colNames=FALSE,overwrite = T)
  
  #SERRF
  setwd(savepath)
  SERRF2(input = "all_data.xlsx")
  
  #SERRF_normalized.csv中去除validate
  SERRF_normalized <- readdata(file =paste0(savepath,"/SERRF_normalized.csv"))
  SERRF_normalized<-rename(SERRF_normalized,c("X"="name"))
  SERRF_normalized<-select(SERRF_normalized,-c(grep("validate",colnames(SERRF_normalized))))
  write.csv(SERRF_normalized,file=paste0(afterpath,"/SERRF_normalized.csv"),row.names=F)
  
  #将分组文件转化为list,为PCA class输入文件
  my_class <- split(class, class$sampleType)
  my_class <- lapply(my_class, function(x) x$label)
  #绘制矫正后图片
  svr_data<-SERRF_normalized%>%remove_rownames()%>%column_to_rownames("name")
  svr_pca<-lmbio::mulstatistics_file(datafile = svr_data,classfile =my_class,mode = "PCA")
  pdf(paste0(afterpath,"/svr_pca.pdf"))
  lmbio::map_mulstatistics_scoremap(filename = svr_pca,
                                    imagetype = NA)
  dev.off()
  png(paste0(afterpath,"/svr_pca.png"))
  lmbio::map_mulstatistics_scoremap(filename = svr_pca,
                                    imagetype = NA)
  dev.off()
  #lmbio::ggplotsave(pp,savepath = beforepath,mapname = "svr_pca")
  Samplecorrplot2(input=svr_data,saminfo=class,savepath=afterpath)
  
  norm_qc <- SERRF_normalized[c(grep("QC",colnames(SERRF_normalized)))]
  qc_rsd<- (apply(norm_qc,1,sd)/apply(norm_qc,1,mean))*100
  norm_sample <- select(SERRF_normalized,-c("name",grep("QC",colnames(SERRF_normalized))))
  sample_rsd<- (apply(norm_sample,1,sd)/apply(norm_sample,1,mean))*100
  rsddata_after <- data.frame(#name =SERRF_normalized[,"name"],
    SERRF_normalized,
    qc_rsd = qc_rsd,
    sample_rsd = sample_rsd
  )	
  write.csv(rsddata_after,paste0(afterpath,"/SERRF_rsd.csv"),row.names=F)
  frequencydraw(frompath = afterpath,file = "SERRF_rsd.csv",savepath=afterpath)
  
  #绘制矫正前图片
  raw_matrix<-select(raw_matrix,-c("No"))%>%remove_rownames()%>%column_to_rownames("label")
  raw_pca <- lmbio::mulstatistics_file(datafile = raw_matrix,
                                       classfile = my_class,
                                       mode = "PCA")
  pdf(paste0(beforepath,"/raw_pca.pdf"))
  lmbio::map_mulstatistics_scoremap(filename = raw_pca,
                                    imagetype = NA)
  dev.off()
  png(paste0(beforepath,"/raw_pca.png"))
  lmbio::map_mulstatistics_scoremap(filename = raw_pca,
                                    imagetype = NA)
  dev.off()
  Samplecorrplot2(input=raw_matrix,saminfo=class,savepath=beforepath)
  
  raw_df<-data.frame(rawdata[,samplestart:ncol(rawdata)])
  rawqc<-raw_df[c(grep("QC",colnames(raw_df)))]
  qc_rsd_raw<- (apply(rawqc,1,sd)/apply(rawqc,1,mean))*100
  rawsample<-select(raw_df,-c(grep("QC",colnames(raw_df))))
  sample_rsd_raw <- (apply(rawsample,1,sd)/apply(rawsample,1,mean))*100
  rsddata_before <- data.frame(name=rawdata[,"name"],
                               qc_rsd = qc_rsd_raw,
                               sample_rsd = sample_rsd_raw,
                               raw_df )
  write.csv(rsddata_before,paste0(beforepath,"raw_rsd.csv"),row.names=F)
  frequencydraw(frompath = beforepath,file = "raw_rsd.csv",savepath=beforepath)
  system("chmod -R 777 ./* ")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio")) 
  
  parser <- ArgumentParser()
  parser$add_argument("-i","--input",default = "./数据矩阵.xlsx", help = "数据矩阵")
  parser$add_argument("-si","--saminfo",default = "./Meta.csv", help = "样本信息数据,包括样本名称、批次、分类和上机顺序")
  parser$add_argument("-sp","--savepath",default = "批次矫正", help = "保存数据路径")
  
  args <- parser$parse_args()
  result <- do.call(what =SERRF_run ,args = args) 
} 
