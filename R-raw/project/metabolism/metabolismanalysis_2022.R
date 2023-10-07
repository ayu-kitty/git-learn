#! /opt/conda/bin/Rscript

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-n","--name",default = "数据矩阵.xlsx", help = "文件名")
  parser$add_argument("-t","--type",default = "非靶向代谢-外来数据", help = "数据来源")
  parser$add_argument("-s","--species",default = "ko", help = "物种")
  parser$add_argument("-c","--ctype",default = "ttest", help = "单变量统计算法")
  parser$add_argument("-fd","--fdr",default = "BH", help = "FDR算法")
  parser$add_argument("-v","--vip",default = 0, type = "double",help = "vip筛选")
  parser$add_argument("-f","--foldchange",default = 0, type = "double",help = "fc筛选")
  parser$add_argument("-pa","--pvalue",default = 0.05, type = "double",help = "p值筛选")
  parser$add_argument("-a","--adjpvalue",default = 0, type = "double",help = "fdr值筛选")
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-z","--zip",default = F, help = "是否压缩",action='store_true')
  
  args <- parser$parse_args()
  
  # print(args)
  
  if(args$type == "非靶向代谢-外来数据"){
    log10L <- F
  }else if(args$type == "非靶向代谢-LCMS"){
    log10L <- F
  }else if(args$type == "非靶向代谢-GCMS"){
    log10L <- T
  }else{
    log10L <- F
  }

  data <- meta::ArrangeInfoFrom(name = args$name,
                                datafrom = "数据矩阵",
                                comparefrom = "比较",
                                saveregistration = F,
                                missvalue = NA,
                                rsd = NA,
                                zeroprocess = NA,
                                type = args$type,
                                species = args$species,
                                log10L = log10L)

  # data <- allstatistics2(data)
  # data <- allstatistics3(data)
  # data <- datastatistics(data)
  # data <- automulstatistics(data)
  data <- run_mulstatisticsana(data)
  data <- meta::getvip(data,statistics=T)
  data <- meta::countdif(data,
                         type = args$ctype,
                         adjust = args$fdr)

  data$info$basic$项目类型 <- "非靶向代谢-外来数据"
  data$info$basic$处理类别 <- "无QC有分析"

  # logfile <- file("log.txt", open = "a")
  # sink(logfile)
  # sink(logfile, type = "message")

  data <- meta::createlist(data,
                           VIP = ifelse(args$vip==0,NA,args$vip),
                           FC = ifelse(args$foldchange==0,NA,args$foldchange),
                           Pvalue = ifelse(args$pvalue==0,NA,args$pvalue),
                           adjPvalue = ifelse(args$adjpvalue==0,NA,args$adjpvalue),
                           family = args$family)

  # 生成报告
  mkreport(savepath = "项目报告",
           type = data$info$basic$项目类型)

  # sink()
  # sink(type = "message")

  if(args$zip){
    zip::zip(zipfile = "项目报告.zip",files = "项目报告")
    unlink(list.files(pattern = "[^z][^i][^p]$",include.dirs = T,all.files = T), recursive = T)
  }
  
}
