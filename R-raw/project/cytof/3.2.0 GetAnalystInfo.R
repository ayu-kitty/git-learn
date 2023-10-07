#!/opt/conda/bin/Rscript

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-ad","--analysis_id", help = "分析编号",required = T)
  parser$add_argument("-sp","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-fn","--filename",default = "项目登记单.xlsx",help = "保存文件名")
  parser$add_argument("-ow","--overwrite",default = F,help = "是否覆盖原始分析单",action='store_true')
  
  args <- parser$parse_args()
  
  result <- do.call(what = GetAnalystInfo,args = args) 
  
}
