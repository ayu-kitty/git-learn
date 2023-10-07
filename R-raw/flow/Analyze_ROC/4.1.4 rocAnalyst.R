#!/opt/conda/bin/Rscript

#' 调用roc分析脚本
#'
#' @param filename 文件名
#' @param mapname 保存图片名
#' @param paramtype roc分析类型：roc、logroc、logrocnum
#' @param saveRoc 逻辑，是否保存分析结果，默认保存
#' @param savelogRoc 逻辑，是否保存分析结果，默认保存
#' @param savelogRocnum 逻辑，是否保存分析结果，默认保存
#' @param number 多roc叠加时特征的个数
#' @param num 随机组合进行逻辑回归roc计算时组合特征的个数
#' @param all 是否保存roc结果数据
#' @param ... 
#'
#' @export
rocAnalyst <- function(filename ,
                       mapname = NULL,
                       paramtype = "roc",
                       number = 1,
                       num = 2,
                       all = T,
                       ...){
  
  getdata <- readdata(filename = filename, 
                      row.names = 1)
  
  if(paramtype == "logrocnum"){
    auto_logrocnum(data = getdata,
                   num = num,
                   all = all,
                   ...)
  }else if(paramtype == "logroc"){
    auto_logroc(data = getdata,
                ...) 
  }else if(paramtype == "roc"){
    auto_roc(data = getdata,
             number = number,
             ...) 
    
  }else{print("输入的绘制类型不支持")}
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",default = "表达数据矩阵", help = "表达数据矩阵文件,如果是xlsx文件以`数据矩阵.xlsx,数据矩阵`形式传参")
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf","html"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 特有参数
  parser$add_argument("-pt","--paramtype", default = "roc", type = "character", help="绘制图片类型：roc、logroc、logrocnum")
  parser$add_argument("-nsr","--nosaveroc", default = T, action = "store_false", 
                      help="是否保存roc结果数据",dest = "saveroc")
  parser$add_argument("-nu","--number", default = 3, type = "integer", help="多roc叠加，要与输入的数据量对应")
  parser$add_argument("-a","--AUC", default = "T", help="是否显示AUC")
  parser$add_argument("-r","--right",default = "T",help = "AUC值是否自动换算成>0.5")
  parser$add_argument("-n","--num", default = 2, type = "integer", help="随机组合进行逻辑回归roc计算")
  parser$add_argument("-al","--all", default = T, action = "store_false", help="是否一起绘制随机组合进行逻辑回归的roc")
  parser$add_argument("-z","--zip",default = F, help = "是否压缩",action='store_true')
  parser$add_argument("-s","--savepath",default = "分析结果", help = "结果输出路径")
  
  args <- parser$parse_args()
  
  zip <- args$zip
  args$zip <- NULL
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  if(args$AUC == "T"){args$AUC <- T}else{args$AUC <- F}
  if(args$right == "T"){args$right <- T}else{args$right <- F}
  
  rocresult <- do.call(what = rocAnalyst, args = args)
  
  if(zip){
    zip::zip(zipfile = "分析结果.zip",files = args$saveptah)
    unlink(list.files(pattern = "[^z][^i][^p]$",include.dirs = T,all.files = T), recursive = T)
  }

}
