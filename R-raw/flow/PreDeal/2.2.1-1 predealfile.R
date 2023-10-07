#!/opt/conda/bin/Rscript

#' 预处理文件
#' 
#' @param filename 保存文件名
#' @param missvarremovepercent 缺失值筛选比例,默认0.5
#' @param missvarremovebygroup 逻辑,是否按组进行缺失值比例计算,使用参数将不按组进行计算
#' @param missvarfillmethod 缺失值填充方式,包含halfmin,valuemin,min,mean,median,knn_var,knn_smp,ppca,bpca,svdlmpute
#' @param rowNorm 归一化方式,包含SumNorm,MedianNorm,QuantileNorm,SamplePQN(ref需提供样本名),GroupPQN(ref需提供组名),CompNorm(ref需提供蛋白或者代谢名),SpecNorm(ref需提供样本名)
#' @param transNorm 转化方式,包含LogNorm,SrNorm,CrNorm
#' @param scaleNorm 标准化方式,包含MeanCenter,AutoNorm,ParetoNorm,RangeNorm
#' @param ref 归一化方式提供参考样本或者特征的参数
#' @param filter 数量比例筛选方式，包含rsd,nrsd,mean,sd,iqr,mad,median,none
#' @param remainnum 特征筛选上限,默认为100000,如需按数量自动输入NULL
#' @param qcFilter 针对QC样本进行rsd筛选
#' @param rsd rsd筛选范围,默认0.3
#'
#' @export
predealfile <- function(filename,
                        missvarremovepercent = 0.5,
                        missvarremovebygroup = T,
                        missvarfillmethod = "halfmin",
                        filter = "rsd",
                        remainnum = 100000,
                        qcFilter = T, 
                        rsd = 30,
                        rowNorm = "NULL",
                        transNorm = "NULL",
                        scaleNorm = "NULL",
                        ref = NULL){
  
  args <- getFunc_Paras()
  
  missvarremovebygroup <- as.logical(missvarremovebygroup)
  if(is.na(missvarremovebygroup)){missvarremovebygroup <- T}
  qcFilter <- as.logical(qcFilter)
  if(is.na(qcFilter)){qcFilter <- T}
  
  if(remainnum == 0){
    remainnum <- NULL
  }
  
  suppressMessages(library("MetaboAnalystR"))
  unlink(x = list.files(pattern = "\\.qs$"))
  mSet <- InitDataObjects("conc", "stat", FALSE)
  mSet <- Read.TextData(mSet,filename, "rowu", "disc")
  mSet <- SanityCheckData(mSet)
  
  mSet <- RemoveMissingPercent2(mSetObj = mSet, 
                                percent = missvarremovepercent,
                                missvarremovebygroup = missvarremovebygroup)
  
  # if(missvarremovepercent > 0){
  #   
  #   args$filename <-  missvaluedel(filename = filename,
  #                                  missvarremovepercent = missvarremovepercent,
  #                                  missvarremovebygroup = missvarremovebygroup)
  #   
  #   args$missvarremovepercent <- 0
  #   do.call(predealfile,args)
  #   return()
  # }else{
  #   args$missvarremovepercent <- 0
  # }
  
  # if(is.null(missvarfillmethod) | missvarfillmethod == "none"){
  #   
  # }else if(missvarfillmethod %in% c("halfmin","valuemin")){
  #   
  #   args$filename <- missvalueprocess(filename = filename,
  #                                     method = missvarfillmethod)
  #   
  #   args$missvarfillmethod <- "none"
  #   do.call(predealfile,args)
  #   return()
  # }
  
  mSet <- ImputeMissingVar2(mSetObj = mSet,
                            method = missvarfillmethod,
                            percent = missvarremovepercent)
  
  mSet <- FilterVariable2(mSetObj = mSet, 
                          filter = filter,
                          remain.num = remainnum,
                          qcFilter = qcFilter, 
                          rsd = rsd)
  
  mSet <- PreparePrenormData(mSet)
  mSet <- Normalization2(mSetObj = mSet,
                         rowNorm = rowNorm,
                         transNorm = transNorm,
                         scaleNorm = scaleNorm,
                         ref = ref)
  
  mSet <- SaveTransformedData(mSet)
  
}

#' @export
getFunc_Paras <- function() {
  pf <- parent.frame()    
  args_names <- ls(envir = pf, all.names = TRUE, sorted = FALSE)
  if("..." %in% args_names) {
    dots <- eval(quote(list(...)), envir = pf)
  }  else {
    dots = list()
  }
  args_names <- sapply(setdiff(args_names, "..."), as.name)
  if(length(args_names)) {
    not_dots <- lapply(args_names, eval, envir = pf) 
  } else {
    not_dots <- list()
  }
  out <- c(not_dots, dots)
  out[names(out) != ""]         
}

