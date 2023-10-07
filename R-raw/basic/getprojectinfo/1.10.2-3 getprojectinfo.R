#!/opt/conda/bin/Rscript

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser1 <- ArgumentParser(Parsername = "parser1",add_help = F)
  group1 <- parser1$add_argument_group('线上登记单获取参数')
  group1$add_argument("-ad","--analysis_id",default = basename(getwd()), help = "分析编号")
  group1$add_argument("-sp","--savepath",default = "./", help = "保存路径")
  # group1$add_argument("-fn","--filename",default = "分析确认单.xlsx",help = "保存文件名")
  group1$add_argument("-ow","--overwrite",default = F,help = "是否覆盖原始分析单",action='store_true')
  
  parser2 <- parser1$newparser(Parsername = "parser2",add_help = F)
  group2 <- parser2$add_argument_group('登记单整理参数')
  # group2$add_argument("-xn","--xlsxname",default = "内部分析单.xlsx",help = "保存文件名")
  group2$add_argument("-tp","--type",default = "",help = "项目类型")
  group2$add_argument("-ks","--keggspecies",default = "ko", help = "物种",required = T)
  
  parser <- parser1$newparser(parents = list("[parser1,parser2]"),Parsername = "parser")
  
  args <- parser$parse_args()
  args1 <- parser1$parse_known_args()[[1]]
  args2 <- parser2$parse_known_args()[[1]]
  
  if(args1$savepath == ""){
    args1$savepath <- NULL
    args2$savepath <- args1$analysis_id
  }else{
    args2$savepath <- args1$savepath
  }
  if(args2$type == ""){args2$type <- NULL}
  # args2$fromname <- args1$filename
  args2$overwrite <- args1$overwrite
  
  # print(args1)
  # print(args2)
  
  result <- do.call(what = GetAnalystInfo,args = args1) 
  result <- do.call(what = ArrangeInfoinpath,args = args2) 
}
