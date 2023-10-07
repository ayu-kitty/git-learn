#!/opt/conda/bin/Rscript

#' 出具oebio体系报告
#'
#' @param oebio python脚本
#' @param config 默认yaml后缀的脚本，支持xlsx为后缀的文件
#' @param logo 鹿明或欧易的logo，默认欧易，修改时填入："lm"
#' @param workd 数据保存的工作目录
#' @param reportdir 报告目录，如需使用自己的脚本，则填yaml所在目录
#' @param type 项目类型，例如：靶向-氨基酸
#' @param date_v 时间版本
#'
#' # @export
# writereport <- function(oebio = "Report.py",
#                         config = "config.yaml",
#                         reportdir = NULL,
#                         type = NULL,
#                         date_v = NULL,
#                         logo ="oe",
#                         workd = "."){
# 
#   suppressMessages(library("tools"))
#   # yaml格式的config文件
#   file_type = file_ext(config) 
#   # 提供脚本的绝对目录，使用模板的python脚本
#   if(is.null(reportdir)){
#     reportdir = "/data/hstore2/home/guj/Cloud/常规报告"
#     
#     if(is.null(date_v)){
#       oebio = dir(reportdir)[max(grep(type, dir(reportdir)))]
#       oebio = file.path(reportdir, file.path(oebio, paste0(unlist(strsplit(type,'-'))[-1], '.py')))
#     }else{
#       oebio = dir(reportdir)[grep(paste(type, date_v, sep = '-', collapse = ''), dir(reportdir))]
#       oebio = file.path(reportdir, file.path(oebio, paste0(unlist(strsplit(type,'-'))[-1], '.py')))
#     }
#     
#     
#     if(file.exists(oebio)){
#       if(!dir.exists(file.path(workd,'src'))){
#         dir.create(file.path(workd,'src'))
#       }
#       
#       if(length(grep("^/",oebio)) == 1){
#         srcpath_company = file.path(dirname(oebio),"../src")
#         srcpath_project = file.path(dirname(oebio),"src")
#         
#         srccopy(srcpath = srcpath_company, dstpath = file.path(workd, "src"))
#         srccopy(srcpath = srcpath_project, dstpath = file.path(workd, "src"))
#     }
#   }
#   
#     # yaml格式的config文件
#     if(file_type=="yaml"){
#       config = paste0(unlist(strsplit(oebio, "py")), "yaml")
#       system(paste("/opt/conda/bin/python3", oebio, "-c", config, "-lg", logo, "-d", workd, sep=" "), intern = T)
#     }else if(file_type=="xlsx"){
#       # 其他用于config的文件需要固定在python脚本中
#       srccopy(srcpath = paste0(unlist(strsplit(oebio, "py")), file_type), dstpath = file.path(workd))
#       system(paste("/opt/conda/bin/python3", oebio, sep=" "), intern = T)
#     }
#   }else{
#     
#     # 不使用模板的python脚本，则报告所需的文件需固定在python脚本中
#     system(paste("/opt/conda/bin/python3", oebio, sep=" "), intern = T)
#   }
#   
# }

#' 生成报告
#'
#' @param savepath 保存路径
#' @param type 报告类型
#' @param ... 
#' @param mode 公司模板
#'
#' @export
mkreport <- function(savepath = "./",
                     type = "",
                     mode = "LM",
                     force = T,
                     CLOUD=F,
                     ...) { 
  
  wd <- getwd()
  
  setwddir(filename = savepath,force = F)
  
  print("生成报告中")
  reportpath <- selectfile(path = "report",file = type,
                           ...)
  
  copydir(from = reportpath,to = "./")
  # srcpath <- paste0(dirname(reportpath),"/src")
  # if(file.exists(srcpath)){
  #   copydir(from =srcpath,to = "src")
  # }
  base::print(force)
  system(paste("/opt/conda/bin/python report.py ", '-m' ,mode,'-cl',ifelse(CLOUD,"T","F"), sep=" "))
  
  # runpython(path = "report.py")
  if(force){
    rmfile <- list.files(path = "./",pattern = ".yaml$|.py$")
    unlink(x = rmfile)
    unlink(x = "src",recursive = T,force = T)
  }
  
  setwd(wd)
}

#' 出具oebio体系报告
#'
#' @param logo 鹿明或欧易的logo，默认欧易，修改时填入："lm"
#' @param type 项目类型，例如：靶向-氨基酸
#' 
#' @export
writereport <- function(type = NULL,logo = "LM"){
  
  mode <- toupper(logo)
  
  mkreport(type = type,mode = mode)
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-s","--savepath",help = "报告路径",required = T)
  parser$add_argument("-t","--type",help = "报告模板",required = T)
  parser$add_argument("-m","--mode",default = "LM", help = "公司模板")
  parser$add_argument("-sm","--selectmode",default = "!", help = "模板选择模式,可填写!(重要性),#(测试),''(时间),版本")
  parser$add_argument("-v","--version",default = lmbio::lmbioversion, help = "大版本填写")
  parser$add_argument("-nf","--nforce",default = T, action = "store_false",help = "是否强制删除模板",dest = "force")
  
  args <- parser$parse_args()
  
  base::print(args$force)
  
  result <- do.call(what = mkreport,args = args) 
  
}
