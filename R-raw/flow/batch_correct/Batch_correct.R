#!/opt/conda/bin/Rscript
#' @export
Batch_correct <- function(input = "表达矩阵.xlsx",
                          saminfo = "批次信息.xlsx",
                          savepath = "批次矫正",
                          method = "statTarget",
                          minqc = 0,
                          minsample = 0,
                          Frule = 0.8, 
                          MLmethod = "QCRFSC", 
                          QCspan = 0,
                          imputeM = "KNN",
                          ...){
  
  wd <- getwd()
  
  if(method=="Combat"){
    result <- ComBat_batch(input = input,saminfo = saminfo,savepath = savepath,...)
  }else if(method=="SERRF"){
    result <- SERRF_run(input = input,saminfo = saminfo,savepath = savepath,...)
  }else if(method=="MetNormalizer"){
    result <- MetNormalizer(input = input,saminfo = saminfo,savepath = savepath,minqc = minqc,minsample = minsample,...)
  }else if(method=="statTarget"){
    result <- Stattarget(input = input,saminfo = saminfo,savepath = savepath,Frule = Frule,MLmethod = MLmethod,QCspan = QCspan,imputeM = imputeM,...)
  }else{
    return(input)
  }
  
  setwd(wd)
  return(result)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-i","--input",default = "./数据矩阵.xlsx", help = "数据矩阵")
  parser$add_argument("-sp","--savepath",default = "批次矫正", help = "保存数据路径")
  parser$add_argument("-si","--saminfo",default ="批次信息.xlsx", help = "包含代谢信息的文件，包括样本名称、批次、分类和上机顺序")
  parser$add_argument("-m","--method",default = "statTarget", help = "去批次的方法：Combat(针对蛋白数据),SERRF,MetNormalizer,statTarget")
  parser$add_argument("-mc","--minqc",default = 0, help = "MetNormalizer参数：最少qc样本数")
  parser$add_argument("-ms","--minsample",default = 0, help = "MetNormalizer参数：最少sample样本数")
  parser$add_argument("-fr","--Frule",default = 0.8, help = "80%过滤原则")
  parser$add_argument("-ml","--MLmethod",default = "QCRFSC", help = "基于QC的信号校正的机器学习方法,QCRFSC:随机森林,QCRLSC:LOESS信号")
  parser$add_argument("-qs","--QCspan",default = 0, help = "QCRLSC的平滑参数，用于控制QCRLSS方法中的偏差-方差平衡,为0表示广义交叉验证,防止过拟合")
  parser$add_argument("-im","--imputeM",default = "KNN", help ="缺失值填充方法")
  
  args <- parser$parse_args()
  result <- do.call(what = Batch_correct,args = args) 
}
