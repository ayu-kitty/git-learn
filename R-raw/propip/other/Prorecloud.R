#!/opt/conda/bin/Rscript

#' 云平台出报告函数
#'
#' @param rawpath 
#' @param reportpath 
#' @param mode 
#' @param cloud 
#' @param ... 
#' @export
Prorecloud<-function(rawpath="./",reportpath="差异分析结果/",cloud="F",mode="LM",...){
  library(stringr,dplyr)
  difdata <- readxlsx(paste0(reportpath,"/result/2.Different_expression/差异表达矩阵.xlsx"), sheet = 1)
  if(!"GO_term" %in% names(difdata)){
    diffanno(reportpath = reportpath)
  }
  system("cp -r /data/hstore3/public/propip/report/云平台/* ./")
  system("chmod -R 777 src")
  # 压缩keggmap
  if(dir.exists("差异分析结果/result/3.Enrichment/KEGG_map/")){
    path0=getwd()
    setwd("差异分析结果/result/3.Enrichment/")
    system("zip -qrm KEGG_map.zip KEGG_map/*")
    setwd(path0)
  }
  # 正式报告
  system(paste0("python Report.py -m ",mode," -c ",cloud," >/dev/null"))
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  parser$add_argument("-m","--mode", default = "LM",help = "报告模板类型，包括LM、OE、QD、SG、HY、GY、YZ、unlogo，默认为LM")
  parser$add_argument("-c","--cloud", default = "F",help = "报告是否云交付，默认为F")
  parser$add_argument("-rd","--rawpath", default = "./",help = "存放出报告文件路径，默认为./")
  parser$add_argument("-rp","--reportpath", default = "差异分析结果/",help = "报告出具路径，默认为差异分析结果/")
  args <- parser$parse_args()
  Prorecloud <- do.call(Prorecloud,args = args)
}
