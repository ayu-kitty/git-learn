#! /opt/conda/bin/Rscript

#' @export 
flow_organize_predealdata<-function(rawfile = "数据矩阵.xlsx",
                                    datafile = "oecloud/rawdata/rawdatafile.txt",
                                    infofile = "oecloud/rawdata/infofile.txt",
                                    classfile = "oecloud/rawdata/classfile.yaml",
                                    classtypefile = "classtype.xlsx",
                                    # comparefile = "oecloud/rawdata/compare.yaml",
                                    processsavepath = "oecloud/predeal",
                                    resultfile = "oecloud/rawdata/datafile.txt",
                                    missvarremovepercent = "0.5",
                                    missvarremovebygroup = T,
                                    missvarfillmethod = "halfmin", 
                                    rowNorm = "NULL", 
                                    transNorm = "NULL",
                                    scaleNorm = "NULL",
                                    filter = "mean",
                                    remainnum = 100000,
                                    qcFilter = T,
                                    rsd = 30){
  organizedata (rawfile = rawfile,
                datafile = datafile,
                infofile = infofile,
                classfile = classfile,
                classtypefile = classtypefile,
                keylist = 1)
  pre_data <- predealdataforfile (datafile = datafile,
                                  classfile = classfile,
                                  processsavepath = processsavepath,
                                  resultfile = resultfile,
                                  missvarremovepercent = missvarremovepercent,
                                  missvarremovebygroup = missvarremovebygroup,
                                  missvarfillmethod = missvarfillmethod, 
                                  rowNorm = rowNorm, 
                                  transNorm = transNorm,
                                  scaleNorm = scaleNorm,
                                  filter = filter,
                                  remainnum = remainnum,
                                  qcFilter = qcFilter,
                                  rsd = rsd)
  infodata <- readdata(infofile)
  #prefile <-merge(pre_data,infodata,by = "FeatureID",all.x = T)
  prefile <- cbind(infodata,pre_data) %>% select(-FeatureID)
  savetxt(data = prefile,filename = "./predatafile.txt",row.names = F,quote = F)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  suppressMessages(library(tidyr))
  suppressMessages(library(dplyr))
  parser <- ArgumentParser()
  
  parser$add_argument("-rf","--rawfile",default = "数据矩阵.xlsx", help = "数据矩阵")
  parser$add_argument("-s","--resultfile",default = "./oecloud/rawdata/datafile.txt", help = "结果存储路径")
  parser$add_argument("-mrp","--missvarremovepercent",default = 0.5, type= "double",
                      help = "缺失值筛选比例,默认0.5")
  parser$add_argument("-mrg","-nmissvarremovebygroup",default = "T", 
                      help = "是否按组进行缺失值比例计算,使用参数将不按组进行计算")
  parser$add_argument("-mvf","--missvarfillmethod",default = "halfmin", 
                      help = "缺失值填充方式,包含none,mean_half,na_half,halfmin,valuemin,min,mean,median,knn_var,knn_smp,ppca,bpca,svdlmpute",
                      choices = c("none","mean_half","na_half","halfmin","valuemin","min","mean","median","knn_var","knn_smp","ppca","bpca","svdlmpute"))
  parser$add_argument("-rn","--rowNorm",default = "NULL", 
                      help = "归一化方式,包含NULL,MedianMeanNorm,SumNorm,MedianNorm,QuantileNorm,SamplePQN(ref需提供样本名),GroupPQN(ref需提供组名),CompNorm(ref需提供蛋白或者代谢名),SpecNorm(ref需提供样本名)",
                      choices = c("NULL","SumNorm","MedianMeanNorm","MedianNorm","QuantileNorm","SamplePQN","GroupPQN","CompNorm"))
  parser$add_argument("-tn","--transNorm",default = "NULL", help = "转化方式,包含LogNorm,Log2minNorm,SrNorm,CrNorm",
                      choices=c("NULL","LogNorm","Log2minNorm","SrNorm","CrNorm"))
  parser$add_argument("-sn","--scaleNorm",default = "NULL", help = "标准化方式,包含MeanCenter,AutoNorm,ParetoNorm,RangeNorm",
                      choices = c("NULL","MeanCenter","AutoNorm","ParetoNorm","RangeNorm"))
  
  parser$add_argument("-fi","--filter",default = "mean", 
                      help = "数量筛选方式,默认为mean,包含rsd,nrsd,mean,sd,iqr,mad,median,none",
                      choices = c("rsd","nrsd","mean","sd","iqr","mad","median","none"))
  parser$add_argument("-re","--remainnum",default = 100000, type= "integer",help = "特征筛选上限,默认为100000,输入0为自动筛选数量")
  parser$add_argument("-qf","--qcFilter",default = "T",help = "是否进行QC样本的rsd筛选")
  parser$add_argument("-rs","--rsd",default = 30,type= "integer", help = "QC样本rsd筛选范围,默认30")
  
  args <- parser$parse_args()
  
  args <- args[!unlist(lapply(args,is.null))]
  
  flowargs <- do.call(what = flow_organize_predealdata,args = args)
}