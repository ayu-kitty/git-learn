#!/opt/conda/bin/Rscript

#' @export
all_get_diff_ana_data <- function(rdspath = "rawdata/oecloud/diff_ana_vip/",
                                  rdsname = list.files(path = rdspath,pattern = "\\.rds$",full.names = T),
                                  savefilename = "差异代谢物(未筛选).xlsx",
                                  savedata = T,
                                  datafrom = "diffdata",
                                  ...){
  
  if (file.exists(savefilename)) {
    wb <- openxlsx::loadWorkbook(savefilename)
  } else {
    wb <- openxlsx::createWorkbook()
  }
  
  for ( i in 1:length(rdsname)) {
    wb <- get_diffana_data(data = rdsname[i],
                           wb = wb,
                           savedata = F,
                           datafrom = datafrom,
                           ...)
  }
  
  if(savedata){
    
    savewb(wb = wb,
           filename = savefilename,
           overwrite = TRUE)
    
  }
  
}

#' 获取差异分析的所有数据
#'
#' @param data obj
#'
#' @export
get_diffana_data <- function(data,
                             datafrom = "diffdata",
                             sort = c("p-value","log2FoldChange","VIP"),
                             decreasing = c(F,F,T),
                             savedata = F,
                             savefilename = "差异代谢物(未筛选).xlsx",
                             wb = openxlsx::createWorkbook(),
                             needlist = NULL,
                             adddata = T,
                             ...) {
  
  data <- readdata(data)
  
  diffdata <- data[["result"]][[datafrom]]
  rawdata <- data[["args"]][["data"]]
  infodata <- data[["args"]][["info"]]
  
  sheetname <- gsub(pattern = "/",replacement = "-vs-",x = data[["args"]][["group"]])
  
  diffdata <- merge(infodata,diffdata,by = 0)
  row.names(diffdata) <- diffdata[,1]
  diffdata <- diffdata[,-1,drop = F]
  
  diffdata <- merge(diffdata,rawdata,by = 0)
  row.names(diffdata) <- diffdata[,1]
  diffdata <- diffdata[,-1,drop = F]
  
  if(!is.null(needlist)){
    
    for ( i in 1:length(needlist)) {
      if(needlist[i] %in% colnames(diffdata)){
        diffdata <- diffdata[!is.na(diffdata[,needlist[i]]),]
      }
    }
    
  }
  
  if(!is.null(sort)){
    for ( i in 1:length(sort)) {
      if(sort[i] %in% colnames(diffdata)){
        diffdata <- diffdata[order(diffdata[,sort[i]],decreasing = decreasing[i]),]
      }
    }
  }
  
  if(savedata){
    savexlsx_fc(data = diffdata,
                filename = savefilename,
                sheet = sheetname)
  }else{
    
    if(!is.null(wb)){

      wb <- addsheet(data = diffdata,
                     wb = wb,
                     sheet = sheetname,
                     sheetmoudle = addsheet_fc)
      
      return(wb)
      
    }else{
      diffdata
    }
    
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-rp","--rdspath",default = "oecloud/diff_filter/", help = "rds数据路径")
  parser$add_argument("-r","--rdsname",
                      default = NULL, 
                      help = "中间过程数据位置",
                      nargs = "+")
  parser$add_argument("-s","--savefilename",default = "差异代谢物(未筛选).xlsx", help = "保存结果")
  
  args <- parser$parse_args()
  
  if(is.null(args$rdsname)){
    args$rdsname <- NULL
  }
  
  flowargs <- do.call(what = all_get_diff_ana_data,args = args)
  
}
