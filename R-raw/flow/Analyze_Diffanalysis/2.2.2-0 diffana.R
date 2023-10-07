#!/opt/conda/bin/Rscript

#' @export
all_diff_ana <- function(comparefile = "oecloud/rawdata/compare.yaml",
                         group = NULL,
                         ...) {
  
  compare <- readdata(filename = comparefile)
  if(is.null(compare$params$all_compare)){
    return()
  }
  compare <- compare$params$all_compare
  
  if(!is.null(group)){
    
    compare <- compare[group]
    
  }
  
  for ( i in 1:length(compare)) {
    
    args <- list(...)
    
    args$group <- compare[[i]]$rawname
    
    if(!("paired" %in% names(args))){
      args$paired <- compare[[i]]$paired
    }
    
    if(!("p_adjust_method" %in% names(args))){
      args$`p_adjust_method` <- compare[[i]]$`p_adjust_method`
    }
    
    if(!("p_method" %in% names(args))){
      args$`p_method` <- compare[[i]]$`p_method`
    }
    
    # base::print(args)
    
    result <- do.call(what = diffana,args = args)
    
  }
  
}

#' 差异分析
#'
#' @param class 
#' @param compare 
#' @param data 
#' @param p_method 
#' @param p_adjust_methods 
#' @param paired 
#' @param rdsname 
#' @param mulstatisticspath 
#' @param mulstatisticsrds 
#' @param mulstatisticsdata 
#' @param diffsavepath 
#' @param diffrdsname 
#' @param log 
#' @param ... 
#'
#' @export
diffana <- function(group,
                    datafile = "oecloud/rawdata/datafile.txt",
                    classfile = "oecloud/rawdata/classfile.yaml",
                    infofile = "oecloud/rawdata/infofile.txt",
                    p_method = "t.test",
                    p_adjust_methods = "BH",
                    paired = F,
                    mulstatisticsrds = paste0("./oecloud/mulstatisticsanalyst/",
                                              "mulstatistic_OPLS-DA-",
                                              gsub(pattern = "\\/",
                                                   replacement = "-vs-",
                                                   x = group),
                                              ".rds"),
                    savediffrds = paste0("./oecloud/diff_ana/diff_ana-",
                                         gsub(pattern = "\\/",
                                              replacement = "-vs-",
                                              x = group),
                                         # if(paired){"~paired"}else{NULL},
                                         ".rds"),
                    log = F,
                    ...){
  # cat(group)
  # datafile <- "oecloud/rawdata/datafile.txt"
  # classfile <- "oecloud/rawdata/classfile.yaml"
  # infofile = "oecloud/rawdata/infofile.txt"
  # group <- "Nicotine-High/Nicotine-Low/Saccharin"
  # log = F
  
  data <- readdata(filename = datafile,row.names = 1)
  class <- readdata(filename = classfile)
  info <- readdata(filename = infofile,row.names = 1)
  
  compare <- group
  # base::print(group)
  # 数据处理
  group2 <- unlist(strsplit(compare,split = "\\/"))
  anaclass <- class[group2]
  if(any(duplicated(Reduce(c,anaclass)))){
    stop(paste0(name,"组中有重复样本"))
  }
  singleclass <- data.frame(Group = NULL)
  for ( i in 1:length(anaclass)) {
    newclass <- anaclass[[i]]
    if(!all(newclass %in% colnames(data))){
      warning(paste0("注意:",paste(newclass[!(newclass %in% colnames(data))],collapse = ";"),"不在数据中"),immediate. = T)
      newclass <- newclass[newclass %in% colnames(data)]
    }
    newsingleclass <- data.frame(Group = rep(names(anaclass)[i],length(newclass)),
                                 check.names = F,
                                 stringsAsFactors = F)
    rownames(newsingleclass) <- newclass
    singleclass <- rbind(singleclass,newsingleclass)
  }
  class <- singleclass[,1]
  data <- data[,row.names(singleclass),drop=F]
  
  print("进行差异分析")
  args <- getFunc_Paras()
  
  # starttime <- lubridate::now(tzone = "Asia/Shanghai")
  
  diffdata <- fc_ana(data = data,
                     class = singleclass,
                     compare = compare,
                     log = log) 
  
  diffdata2 <- p_ana(data = data,
                     class = singleclass,
                     paired = paired,
                     p_method = p_method,
                     p_adjust_methods = p_adjust_methods,
                     ...)
  
  if(!is.null(diffdata2)){
    diffdata <- merge(diffdata,diffdata2,by = 0,sort = F,all = T)
    row.names(diffdata) <- diffdata[,1]
    diffdata <- diffdata[,-1]
  }
  
  # cat(mulstatisticsrds)
  if(file.exists(mulstatisticsrds)){
    mulstatisticsdata <- readdata(mulstatisticsrds)
    if("mulstatistics" %in% class(mulstatisticsdata)){
      if(!is.null(mulstatisticsdata[["statistics"]])){
        vipdata <- getvip(mulstatisticsdata)
        vipdata <- vipdata$vipdata[,"VIP",drop = F]
      }else{
        vipdata <- NULL
      }
    }else if("mulstatistics-both" %in% class(mulstatisticsdata)){
      if(!is.null(mulstatisticsdata$lcresult[["statistics"]])){
        vipdatalc <- getvip(data = mulstatisticsdata$lcresult)
        vipdatalc <- vipdatalc$vipdata[,"VIP",drop = F]
      }else{
        vipdatalc <- NULL
      }
      if(!is.null(mulstatisticsdata$gcresult[["statistics"]])){
        vipdatagc <- getvip(data = mulstatisticsdata$gcresult)
        vipdatagc <- vipdatagc$vipdata[,"VIP",drop = F]
      }else{
        vipdatagc <- NULL
      }
      
      if(is.null(vipdatalc) & is.null(vipdatagc)){
        vipdata <- NULL
      }else if(is.null(vipdatalc)){
        vipdata <- vipdatagc
      }else if(is.null(vipdatagc)){
        vipdata <- vipdatalc
      }else{
        vipdata <- rbind(vipdatalc,vipdatagc)
      }
      
    }else{
      vipdata <- NULL
    }
  }else{
    vipdata <- NULL
  }
  
  if(!is.null(vipdata)){
    diffdata <- merge(vipdata,diffdata,by = 0,sort = F,all = T)
    row.names(diffdata) <- diffdata[,1]
    diffdata <- diffdata[,-1]
  }
  
  # endtime <- lubridate::now(tzone = "Asia/Shanghai")
  
  args <- list(path = getwd(),
               # starttime = starttime,
               # endtime = endtime,
               fun = diffana,
               args = args,
               result = list(diffdata = diffdata))
  
  class(args) <- c(class(args),"diff_ana")
  
  saverds(data = args,
          filename = savediffrds)
  
  print(paste0("~",group,"的",p_method,"方法单变量计算完成"))
  
  # return(args)
  return(savediffrds)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-df","--datafile",
                      default = "oecloud/rawdata/datafile.txt",help = "数据矩阵文件",nargs = "+")
  parser$add_argument("-cf","--classfile",
                      default = "oecloud/rawdata/classfile.yaml", help = "分组文件",nargs = "+")
  parser$add_argument("-if","--infofile",
                      default = "oecloud/rawdata/infofile.txt", help = "分组文件",nargs = "+")
  parser$add_argument("-co","--comparefile",
                      default = "oecloud/rawdata/compare.yaml", help = "比较文件",nargs = "+")
  parser$add_argument("-g","--group",default = NULL, help = "不输入比较文件，直接输入比较组",nargs = "+")
  parser$add_argument("-pm","--p_method",default = NULL, 
                      help = "p值计算方法,包含t.test,f.test,fttest,wilcox.test,oneway.test,kruskal.test")
  parser$add_argument("-pa","--p_adjust_methods",default = NULL, 
                      help = "adjp计算方法,包含holm,hochberg,hommel,bonferroni,BH,BY,fdr,none")
  parser$add_argument("-l","--log",default = "F", 
                      help = "数据是否经过log标准化")
  parser$add_argument("-mr","--mulstatisticsrds",default = NULL, 
                      help = "多元统计中间过程文件存储路径")
  parser$add_argument("-sr","--savediffrds",default = NULL, 
                      help = "结果文件存储路径")
  
  args <- parser$parse_args()
  
  if(is.null(args$mulstatisticsrds)){args$mulstatisticsrds <- NULL}
  if(is.null(args$savediffrds)){args$savediffrds <- NULL}
  if(is.null(args$p_method)){args$p_method <- NULL}
  if(is.null(args$p_adjust_methods)){args$p_adjust_methods <- NULL}
  
  args$log <- as.logical(args$log)
  
  flowargs <- do.call(what = all_diff_ana,args = args)
  
}
