#!/opt/conda/bin/Rscript

#' 预处理数据
#' 
#' @param data 数据
#' @param class 分组
#' @param ... 见`predealfile`
#'
#' @return 返回形式 
#' @export
predealdata <- function(data,
                        class = NULL,
                        qcFilter = T,
                        ...){
  
  qcFilter <- as.logical(qcFilter)
  if(is.na(qcFilter)){qcFilter <- T}
  
  if(is.null(class)){
    warning("未给分组信息，自动分组",immediate. = T)
    class <- data.frame(class = rep("1",dim(data)[2]))
    row.names(class) <- colnames(data)
  }
  
  args <- saverawdata(data = data,class = class)
  
  if("QC" %in% class[,1]){
    predealfile(filename = args$filename,
                qcFilter = qcFilter,
                ...) 
  }else{
    predealfile(filename = args$filename,
                qcFilter = F,
                ...)
  }
  
  predata <- readpredata(colname = colnames(data))
  
  return(predata)  
}

#' 根据文件预处理数据
#'
#' @param datafile 数据文件 
#' @param classfile 分组文件
#' @param ... 见`predealfile`
#'
#' @export
predealdataforfile <- function(datafile = "./oecloud/rawdata/rawdatafile.txt",
                               classfile = "./oecloud/rawdata/classfile.yaml",
                               processsavepath = "./oecloud/predeal/",
                               resultfile = "./oecloud/rawdata/datafile.txt",
                               ...){
  options(warn = -1)
  # datafile <- c("数据矩阵.xlsx","数据矩阵")
  # classfile <- c("数据矩阵.xlsx","分组")
  data <- readdata(datafile)
  # rownames(data) <- data[,1]
  
  if(is.null(classfile)){
    stop("未给分组信息")
    # warning("未给分组信息，自动分组",immediate. = T)
    # info <- data[,1,drop = F]
    # data <- data[,-1,drop = F]
    # class <- data.frame(class = rep("1",dim(data)[2]))
    # row.names(class) <- colnames(data)
  }else{
    class <- readdata(classfile,row.names = 1)
    if(!is.data.frame(class)){
      sampledata <- class
      class <- NULL
      
      for ( i in 1:length(sampledata)) {
        
        class2 <- data.frame("Sample" = sampledata[[i]],
                             "Group" = names(sampledata)[i],
                             row.names = sampledata[[i]])
        
        class <- rbind(class,class2)
      }
      class <- class[!duplicated(class[,1]),]
      class <- class[,-1,drop = F]
    }
    info <- data[,!(colnames(data) %in% row.names(class)),drop = F]
    data <- data[,colnames(data) %in% row.names(class),drop = F]
  }
  
  predata <- runinpath(path = processsavepath,
                       moudle = predealdata,
                       moudlename = "数据预处理",
                       data = data,
                       class = class,
                       ...)
  
  predata <- merge(info,predata,by = 0,sort = F)[,-1,drop = F]
  
  if(!is.null(resultfile)){
    if(datafile[1] %in% resultfile){
      file.rename(from = datafile[1],to = paste0(dirname(datafile[1]),"/",
                                                 gsub(pattern = "\\.",replacement = "-raw.",x = basename(datafile[1]))))
    }else if(file.exists(resultfile)){
      file.rename(from = resultfile,to = paste0(dirname(resultfile),"/",
                                                gsub(pattern = "\\.",replacement = "-raw.",x = basename(resultfile))))
    }
    
    savetxt(data = predata,
            filename = resultfile) 
    
  }
  
  return(predata)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))

  parser <- ArgumentParser()
  
  parser$add_argument("-df","--datafile",default = "oecloud/rawdata/rawdatafile.txt", 
                      help = "数据矩阵文件,如果是xlsx文件以`数据矩阵.xlsx 数据矩阵`形式传参",
                      nargs = "+")
  parser$add_argument("-cf","--classfile",default = "./oecloud/rawdata/classfile.yaml", 
                      help = "分组文件,如果是xlsx文件以`数据矩阵.xlsx 分组`形式传参",
                      nargs = "+")
  parser$add_argument("-s","--resultfile",default = "./oecloud/rawdata/datafile.txt", help = "结果存储路径")
  parser$add_argument("-mrp","--missvarremovepercent",default = 0.5, type= "double",
                      help = "缺失值筛选比例,默认0.5")
  parser$add_argument("-mrg","--missvarremovebygroup",default = "T", help = "是否按组进行缺失值比例计算,使用参数将不按组进行计算")
  parser$add_argument("-mvf","--missvarfillmethod",default = "halfmin", 
                      help = "缺失值填充方式,包含none,mean_half,na_half,halfmin,valuemin,min,mean,median,knn_var,knn_smp,ppca,bpca,svdlmpute",
                      choices = c("none","mean_half","na_half","halfmin","valuemin","min","mean","median","knn_var","knn_smp","ppca","bpca","svdlmpute"))
  parser$add_argument("-rn","--rowNorm",default = "NULL", 
                      help = "归一化方式,包含NULL,SumNorm,MedianMeanNorm,MedianNorm,QuantileNorm,SamplePQN(ref需提供样本名),GroupPQN(ref需提供组名),CompNorm(ref需提供蛋白或者代谢名),SpecNorm(ref需提供样本名)",
                      choices = c("NULL","SumNorm","MedianMeanNorm","MedianNorm","QuantileNorm","SamplePQN","GroupPQN","CompNorm"))
  parser$add_argument("-tn","--transNorm",default = "NULL", help = "转化方式,包含LogNorm,Log2minNorm,SrNorm,CrNorm",
                      choices=c("NULL","LogNorm","Log2minNorm","SrNorm","CrNorm"))
  parser$add_argument("-sn","--scaleNorm",default = "NULL", help = "标准化方式,包含MeanCenter,AutoNorm,ParetoNorm,RangeNorm",
                      choices = c("NULL","MeanCenter","AutoNorm","ParetoNorm","RangeNorm"))
  parser$add_argument("-r","--ref",default = "NULL", help = "归一化方式提供参考样本或者特征的参数")
  parser$add_argument("-fi","--filter",default = "mean", 
                      help = "数量筛选方式,默认为mean,包含rsd,nrsd,mean,sd,iqr,mad,median,none",
                      choices = c("rsd","nrsd","mean","sd","iqr","mad","median","none"))
  parser$add_argument("-re","--remainnum",default = 100000, type= "integer",
                      help = "特征筛选上限,默认为100000,输入0为自动筛选数量")
  parser$add_argument("-qf","--qcFilter",default = "T",help = "是否进行QC样本的rsd筛选")
  parser$add_argument("-rs","--rsd",default = 30,type= "integer", help = "QC样本rsd筛选范围,默认30")
  parser$add_argument("-pp","--processsavepath",default = "./oecloud/predeal/", help = "中间过程文件保存位置")
  
  args <- parser$parse_args()
  
  data <- do.call(what = predealdataforfile,args = args)
  
}