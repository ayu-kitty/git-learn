#!/opt/conda/bin/Rscript

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-s","--savepath",help = "报告路径",required = T)
  parser$add_argument("-m","--mode",default = "LM", help = "公司模板")
  args <- parser$parse_args()
  
  mkreport(savepath = args$savepath,
           type = "质谱流式",
           mode = args$mode)
}
