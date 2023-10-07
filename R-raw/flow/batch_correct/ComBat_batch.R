#!/opt/conda/bin/Rscript

#关于去除批次效应的函数有Combat和Combat_seq两个，而Combat主要面对的是带有小数的数据（比方说芯片数据），基于的是贝叶斯原理；而Combat_seq主打RNA-seq的count数据，基于负二项分布回归
#' @export
ComBat_batch <- function(input = "表达矩阵.xlsx",
                         saminfo = "批次信息.xlsx",
                         savepath = "批次矫正"){
  
  suppressWarnings(library(sva))
  suppressWarnings(library(dplyr))
  suppressWarnings(library(tidyverse))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(reshape))
  if (!file.exists(savepath)) {
    dir.create(savepath)
  }
  
  afterpath <- paste0(savepath,"/校正后/")
  beforepath <- paste0(savepath,"/校正前/")

  cdata <- readdata(input)
  # row.names(cdata) <- paste0("_",1:dim(cdata)[1])
  sif <- readdata(saminfo)
  csif <- select(sif,c("sample","batch"))#读取批次信息
  bd <- cdata[,csif[,1]]
  infodata <- cdata[,!(colnames(cdata) %in% csif[,1])]
  group <- unique(csif$batch)
  data_group <- data.frame(FeatureID = row.names(cdata))
  #统计每个批次的非空样本数
  for (i in 1:length(group)) {
    one_group <- bd[,csif[which(csif$batch %in% group[i]),1]]
    one_result <- as.data.frame(apply(one_group,1,function(x) length(na.omit(x))))
    colnames(one_result) <- group[i]
    data_group <- cbind(data_group,one_result)
  }
  
  #保留每个批次至少有两个非空样本的蛋白，rowMeans对逻辑值进行运算，T=1,F=0
  dat <- data_group[rowMeans(data_group[, 2:ncol(data_group)] >1) == 1, ]
  bd <- bd[dat[,1],]
  bd_df <- data.frame(FeatureID = row.names(bd),bd,check.names = F)
  savetxt(data = cdata,filename = paste0(beforepath,"/data.txt"))
  savetxt(data = bd_df,filename = paste0(beforepath,"/数据矩阵-矫正前.txt"))
  
  #将分组文件转化为list,为PCA class输入文件
  class <- select(sif,c("sample","batch","class"))
  my_class <- split(class, class$class)
  my_class <- lapply(my_class, function(x) x$sample)
  group <- unique(sif$class)
  #矫正前PCA和相关性热图
  raw_pca <- mulstatistics_file(datafile = bd,
                                classfile = my_class,
                                group = group,
                                mode = "PCA")

  map_mulstatistics_scoremap(filename = raw_pca,
                             mapname = "raw_pca",
                             savepath = beforepath,
                             imagetype = "png")

  map_common_corrplot2(filename = t(bd),
                       mapname = "raw_corrplot",
                       savepath = beforepath,
                       imagetype = "png")
  
  #去除批次效应
  modcombat <- model.matrix(~1, data = csif)
  batch <- csif$batch
  combat_edata <- ComBat(dat = bd, 
                         batch = batch, 
                         mod = modcombat,
                         par.prior = TRUE,
                         mean.only = F,
                         prior.plots = F)
  
  if(all(bd >= 0)){
    combat_edata[combat_edata <= 0 ] <- min(combat_edata[combat_edata > 0],na.rm = T)
  }
  
  combat_edata_df <- data.frame(FeatureID=row.names(combat_edata),combat_edata,check.names = F)
  savetxt(data = combat_edata_df,filename = paste0(afterpath,"/数据矩阵-矫正后.txt"))
  
  #矫正后PCA和相关性热图
  svr_pca <- mulstatistics_file(datafile = combat_edata,
                                classfile = my_class,
                                group = group,
                                mode = "PCA")
  
  map_mulstatistics_scoremap(filename = svr_pca,
                             mapname = "svr_pca",
                             savepath = afterpath,
                             imagetype = "png")
  
  map_common_corrplot2(filename = t(combat_edata),
                       mapname = "srv_corrplot",
                       savepath = afterpath,
                       imagetype = "png")
  
  result <- merge(infodata,combat_edata,by = 0)
  row.names(result) <- result[,1]
  result <- result[,-1,drop = F]
  
  savetxt(data = result,filename = paste0(afterpath,"/data.txt"))
  system(paste0("chmod -R 777 ",savepath))
  return(result)
  
  #拼接数据
  #ncol(data)-ncol(combat_edata)->aa
  #cbind(data[,1:aa],combat_edata)->newdata
  #write.xlsx(newdata,paste0(savepath,"/去除批次数据矩阵_anno.xlsx"),sheetName="去除批次数据矩阵")
  #计算rsd
  # combat_rsd <- data.frame(Accession=data[,1])
  # for (i in 1:length(group)){ 
  #   subdata <- combat_edata[c(grep(as.character(group[i]),colnames(combat_edata)))]
  #   rsd <- as.data.frame(apply(subdata,1,sd)/apply(subdata,1,mean))*100
  #   colnames(rsd) <- paste0(group[i],"_rsd")
  #   combat_rsd <- data.frame(combat_rsd,rsd)
  # } 
  # all_data <- data.frame(combat_rsd,combat_edata)
  # savetxt(data = all_data,filename = paste0(afterpath,"/去除批次数据_rsd.txt"))
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-i","--input",default = "表达矩阵.xlsx", help = "表达矩阵")
  parser$add_argument("-si","--saminfo",default = "批次信息.xlsx", help = "样本批次信息数据，包括样本名称、批次、分类")
  parser$add_argument("-sp","--savepath",default = "批次矫正", help = "保存数据路径")
  
  args <- parser$parse_args()
  result <- do.call(what = ComBat_batch,args = args) 
} 
