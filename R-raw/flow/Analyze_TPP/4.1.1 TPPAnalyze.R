#!/opt/conda/bin/Rscript

# 使用说明:http://192.168.10.200:8787/help/library/lmbio/doc/TPPAnalyze.html

#' TPPAnalyze热蛋白质组分析
#'
#' @param filepath 输入文件
#' @param namelist 基因或者蛋白的列名
#' @param feature 根据提供特征列表进行筛选
#' @param maxfeaturenum 最大特征数量
#' @param normalize 是否归一化
#' @param resultPath 结果文件路径
#' @param nCores 使用核心数量
#' @param methods 分析方法
#' @param ... 见TPP::TRresults
#'
#' @export
TPPAnalyze <- function(filepath = "测试数据.xlsx",
                       namelist = "Gene Name",
                       feature = NULL,
                       maxfeaturenum = 3000,
                       normalize = F,
                       resultPath = "TPP分析",
                       nCores = 4,
                       methods = "meltcurvefit",
                       pValFilter = list(minR2 = 0.1, maxPlateau = 0.3),
                       zip = T,
                       ...){
  suppressMessages(library("TPP"))
  
  data <- readdata(filename = filepath,sheet = "数据矩阵")
  data <- data[!duplicated(data[,namelist]),]
  
  if(length(feature) > 0){
    data <- data[data[,namelist] %in% feature,]
  }
  
  if(dim(data)[1] > maxfeaturenum){
    warning(paste0("特征超过",maxfeaturenum,",将多余特征进行后续分析"),immediate. = T)
    data <- data[1:maxfeaturenum,]
  }
  
  group <- readdata(filename = filepath,sheet = "分组")
  anagroup <- group[,c("Group","Condition")]
  anagroup <- anagroup[!duplicated(anagroup),]
  config <- NULL
  anadata <- list()
  for (i in 1:dim(anagroup)[1]) {
    config2 <- data.frame(Experiment = anagroup[i,"Group"],
                          Condition = anagroup[i,"Condition"],
                          check.names = F)
    config3 <- group[group$Group == anagroup[i,"Group"] & group$Condition == anagroup[i,"Condition"], ]
    config4 <- data.frame(config3$Temp)
    colnames(config4) <- anagroup[i,"Group"]
    rownames(config4) <- config3$Name
    config4 <- as.data.frame(t(config4))
    config2 <- cbind(config2,config4)
    config <- rbind(config,config2)
    
    anadata2 <- data[,c(namelist,config3$sample)]
    colnames(anadata2) <- c("ID",paste0("Express_",config3$Name))
    anadata2$qupm <- 1:dim(anadata2)[1]
    anadata2 <- list(anadata2)
    names(anadata2) <- anagroup[i,"Group"]
    anadata <- c(anadata,anadata2)
  }
  
  if("比较" %in% getsheetname(filename = filepath)){
    compare <- readdata(filename = filepath,sheet = "比较")
    
    for (i in 1:length(compare)) {
      comparename <- paste0("ComparisonVT",i) 
      config[config$Experiment %in% unlist(strsplit(x = compare[1,1],split = "\\/")),comparename] <- "x"
    }
  }
  
  runtry <- try({  
    TRresults <- analyzeTPPTR(configTable = config,
                              methods = methods,
                              data = anadata,
                              nCores = nCores,
                              normalize =  normalize,
                              resultPath = resultPath,
                              idVar = "ID",
                              fcStr = "Express_",
                              plotCurves = T,
                              pValFilter = pValFilter,
                              ...)
  },
  silent = F)

  resultfile <- paste0(resultPath,"/results_TPP_TR.xlsx")
  if(file.exists(resultfile)){
    resultdata <- readdata(resultfile)
    meltPointname <- colnames(resultdata)[grep(pattern = "^meltPoint_",x = colnames(resultdata))]
    for (i in 1:length(meltPointname)) {
      resultdata <- resultdata[!is.na(resultdata[,meltPointname[i]]),]
    }
    if(dim(resultdata)[1] > 0){
      createdir(filename = paste0(resultPath,"/Melting_Curves_ok"))
      file.copy(from = paste0(resultPath,"/Melting_Curves/meltCurve_",resultdata$meltcurve_plot,".pdf"),
                to = paste0(resultPath,"/Melting_Curves_ok/"))
    }
  }
  
  write(x = "1. dataObj
存放运行过程中的中间数据

2. Melting_Curves
存放所有融化曲线
  
3. Melting_Curves_ok
挑选两组均存在融化值的曲线
  
4. QCplots.pdf
在不同温度下两组间强度分布
  
5. results_TPP_TR.xlsx
所有数据结果",
        file = paste0(resultPath,"/说明.txt"))
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filepath",default = "测试数据.xlsx", help = "输入文件路径",required = T)
  
  # 基本参数
  parser$add_argument("-n","--namelist",default = "Gene Name", help = "Gene或者Protein的列名")
  parser$add_argument("-m","--maxfeaturenum",default = 3000, type= "integer",help = "最大蛋白特征数量")
  parser$add_argument("-r","--resultPath",default = "TPP分析", help = "结果输出路径")
  parser$add_argument("-nm","--normalize",default = F, help = "是否进行标准化",action='store_true')
  parser$add_argument("-mt","--methods",default = "meltcurvefit", help = "TPP运算算法",
                      choices = c("meltcurvefit","splinefit"))
  parser$add_argument("-z","--zip",default = F, help = "是否压缩",action='store_true')
  
  args <- parser$parse_args()
  
  writeinfo()
  result <- do.call(what = TPPAnalyze,args = args)
  writeinfo(endtime = T)
  
  if(args$zip){
    zip::zip(zipfile = "项目报告.zip",files = args$resultPath)
    unlink(list.files(pattern = "[^z][^i][^p]$",include.dirs = T,all.files = T), recursive = T)
  }
  
}
