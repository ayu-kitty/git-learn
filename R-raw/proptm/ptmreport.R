#!/opt/conda/bin/Rscript
#' 出报告函数
#'
#' @param rawpath 
#' @param reportpath 
#' @param mode 
#' @param cloud 
#' @param ... 
#' @export
ptmreport<-function(rawpath="rawdata/",reportpath="report/",cloud="F",mode="auto",ptmfile=NULL,...){
  library(stringr,dplyr)

  system(paste0("cp -r /data/hstore3/public/propip/report/proptm/* ", rawpath))
  if(mode=="auto"){
    ty<-str_extract_all(ptmfile,pattern = 'LM|QD|OE|SG|HY')[[1]][1]
    if(is.na(ty)){
      mode<-"LM"
    }else mode<-ty
  }
  setwd(rawpath)
  # 压缩keggmap
  if(dir.exists("../report/result/4.Enrichment/KEGG_map/")){
    path0=getwd()
    setwd("../report/result/4.Enrichment/")
    system("zip -qrm KEGG_map.zip KEGG_map")
    setwd(path0)
  }
  # 正式报告
  # system(paste0("python Report.py -m ",mode," -c ",cloud," > ../project_report_stderr.log"))
  system(paste0("python Report.py -m ",mode," -c ",cloud," >/dev/null"))
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  parser$add_argument("-m","--mode", default = "auto",help = "报告模板类型，包括LM、OE、QD、SG、HY、GY、YZ、unlogo，默认为自动判断")
  parser$add_argument("-c","--cloud", default = "F",help = "报告是否云交付，默认为F")
  parser$add_argument("-rd","--rawpath", default = "rawdata/",help = "存放出报告文件路径，默认为rawdata/")
  parser$add_argument("-rp","--reportpath", default = "report/",help = "报告出具路径，默认为./report/")
  parser$add_argument("-pt","--ptmfile", default = "auto",help = "项目编号路径")
  args <- parser$parse_args()
  ptmreport <- do.call(ptmreport,args = args)
}