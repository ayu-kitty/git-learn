#!/opt/conda/bin/Rscript

#' @export
all_get_diff_ana_data_2 <- function(rdspath = "oecloud/diff_filter/",
                                    rdsname = list.files(path = rdspath,pattern = "\\.rds$",full.names = T),
                                    savepath = "./",
                                    datafrom = "filterdata",
                                    savedata = T,
                                    ...){
  
  for ( i in 1:length(rdsname)) {
    
    data <- readdata(rdsname[i])
    
    if(datafrom == "filterdata"){
      filterparam <- paste0(ifelse(data[["result"]][["filterargs"]]$vipfilter == 0,"",paste0("_vip-",data[["result"]][["filterargs"]]$vipfilter)),
                            ifelse(data[["result"]][["filterargs"]]$adjpfilter == 0,
                                   ifelse(data[["result"]][["filterargs"]]$pfilter == 0,"",
                                          paste0("_p-val-",data[["result"]][["filterargs"]]$pfilter)),
                                   paste0("_q-val-",data[["result"]][["filterargs"]]$adjpfilter)),
                            ifelse(data[["result"]][["filterargs"]]$fcfilter == 0,"",paste0("_fc-",data[["result"]][["filterargs"]]$fcfilter)))
    }else{
      filterparam <- ""
    }

    
    savefilename <- paste0(savepath,"/diff-data-",gsub(pattern = "/",replacement = "-vs-",x = data[["args"]][["compare"]]),filterparam,".xls")
    
    wb <- get_diffana_data_2(data = data,
                             datafrom = datafrom,
                             savefilename = savefilename,
                             savedata = savedata,
                             ...)
  }
  
  return(wb)
}

#' 获取差异分析的所有数据
#'
#' @param data obj
#'
#' @export
get_diffana_data_2 <- function(data,
                               datafrom = "filterdata",
                               sort = c("p-value","log2FoldChange","VIP"),
                               decreasing = c(F,F,T),
                               savedata = F,
                               savefilename = "差异代谢物(未筛选).xls",
                               needlist = NULL,
                               adddata = T,
                               ...) {
  
  data <- readdata(data)
  
  diffdata <- data[["result"]][[datafrom]]
  rawdata <- data[["args"]][["data"]]
  infodata <- data[["args"]][["info"]]
  
  sheetname <- gsub(pattern = "/",replacement = "-vs-",x = data[["args"]][["group"]])
  
  diffdata <- merge(infodata,diffdata,by = 0,all = F)
  row.names(diffdata) <- diffdata[,1]
  diffdata <- diffdata[,-1,drop = F]
  
  if(adddata){
    diffdata <- merge(diffdata,rawdata,by = 0,all = F)
    row.names(diffdata) <- diffdata[,1]
    diffdata <- diffdata[,-1,drop = F]
  }

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
    savexls(data = diffdata,
            filename = savefilename)
    return(NULL)
  }else{
    return(diffdata)
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
  parser$add_argument("-s","--savepath",default = "./", help = "保存结果路径")
  parser$add_argument("-d","--datafrom",default = "filterdata", help = "数据来源")
  
  args <- parser$parse_args()
  
  if(is.null(args$rdsname)){
    args$rdsname <- NULL
  }
  
  flowargs <- do.call(what = all_get_diff_ana_data_2,args = args)
  
}
