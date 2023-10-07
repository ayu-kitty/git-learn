
#' @export
getpredealparams <- function(omic = "M",
                             MStype = "auto",
                             name = "rowNorm",
                             value = "auto"){
  
  if(omic %in% c("M","ML")){
    
    if(MStype %in% c("LC")){
      argslist <- list("missvarremovepercent" = 0.5,
                       "missvarremovebygroup" = "T",
                       "missvarfillmethod" = "halfmin",
                       "rowNorm" = "NULL",
                       "transNorm" = "NULL",
                       "scaleNorm" = "NULL",
                       "ref" = "NULL",
                       "filter" = "mean",
                       "remainnum" = 100000,
                       "qcFilter" = "T",
                       "rsd" = 30,
                       "log10L" = "F")
    }else if(MStype %in% c("GC")){
      argslist <- list("missvarremovepercent" = 0.5,
                       "missvarremovebygroup" = "T",
                       "missvarfillmethod" = "halfmin",
                       "rowNorm" = "NULL",
                       "transNorm" = "NULL",
                       "scaleNorm" = "NULL",
                       "ref" = "NULL",
                       "filter" = "mean",
                       "remainnum" = 100000,
                       "qcFilter" = "F",
                       "rsd" = 30,
                       "log10L" = "T")
    }else{
      argslist <- list("missvarremovepercent" = 0.5,
                       "missvarremovebygroup" = "T",
                       "missvarfillmethod" = "halfmin",
                       "rowNorm" = "NULL",
                       "transNorm" = "NULL",
                       "scaleNorm" = "NULL",
                       "ref" = "NULL",
                       "filter" = "mean",
                       "remainnum" = 100000,
                       "qcFilter" = "T",
                       "rsd" = 30,
                       "log10L" = "F")
    }

  }else if(omic %in% c("P","PF")){
    
    if(MStype %in% c("PD","PDM","PT")){
      
      argslist <- list("missvarremovepercent" = 0,
                       "missvarremovebygroup" = "F",
                       "missvarfillmethod" = "none",
                       "rowNorm" = "MedianNorm",
                       "transNorm" = "Log2minNorm",
                       "scaleNorm" = "NULL",
                       "ref" = "NULL",
                       "filter" = "mean",
                       "remainnum" = 100000,
                       "qcFilter" = "F",
                       "rsd" = 30,
                       "log10L" = "F")
      
    }else if(MStype %in% c("LF","PL","SP","PU","auto")){
      
      argslist <- list("missvarremovepercent" = 0.5,
                       "missvarremovebygroup" = "T",
                       "missvarfillmethod" = "mean_half",
                       "rowNorm" = "MedianNorm",
                       "transNorm" = "Log2minNorm",
                       "scaleNorm" = "NULL",
                       "ref" = "NULL",
                       "filter" = "mean",
                       "remainnum" = 100000,
                       "qcFilter" = "F",
                       "rsd" = 30,
                       "log10L" = "F")
      
    }else{
      stop(paste0("未找到",MStype,"项目子类型"))
    }
    
  }else{
    stop(paste0("未找到",omic,"项目类型"))
  }
  
  if(value == "auto"){
    value <- argslist[[name]]
  }else if(value < 0){
    value <- argslist[[name]]
  }
  
  return(value)
}

