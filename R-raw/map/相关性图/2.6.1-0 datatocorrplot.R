#!/opt/conda/bin/Rscript

#' 计算相关性
#' 
#' @param ci 是否计算置信区间
#' @param xname x数据名称
#' @param yname y数据名称
#' @param corfilter 相关性筛选标准
#' @param corfiltertype 相关性筛选方法,包含+-,+,-
#' @param pfilter 显著性筛选标准
#' @param adjust 显著性校正算法,包括BH,holm,hochberg,hommel,bonferroni,BY,fdr,none
#' @param transx 是否转置x数据
#' @param transy 是否转置y数据
#' @param use 相关性计算方法,包括pairwise,complete
#' @param ... 见psych::corr.test
#' @param filename 
#' @param data 
#' @param filenamey 
#' @param cormethod 
#' @param mapmoudle 
#' @param mapname 
#' @param savepath 
#'
#' @export 
datatocorrplot <- function(filename = "1",
                           data = readdata(filename = filename,row.names = 1),
                           filenamey = NULL,
                           ci = FALSE,
                           xname = "Featurex",
                           yname = "Featurey",
                           corfilter = 0.95,
                           corfiltertype = "+-",
                           pfilter = 0.05,
                           adjust = "none",
                           transx = T,
                           transy = T,
                           cormethod = "pearson",
                           use = "pairwise",
                           mapmoudle = corrplot,
                           mapname = "corrplot",
                           savepath = "./",
                           savecorr = T,
                           ...){
  
  x <- data
  y <- readdata(filename = filenamey,row.names = 1)
  
  if(is.null(data)){
    
    savetxt(data = "数据为空,不进行相关性分析",
            filename = paste0(savepath,"/说明.txt"))
    
    return()
    
  }
  
  if(dim(data)[1] < 3 & is.null(y)){
    
    savetxt(data = "数据太少,不进行相关性分析",
            filename = paste0(savepath,"/说明.txt"))
    
    return()
    
  }
  
  if(dim(data)[2] < 3){
    
    savetxt(data = "数据太少,不进行相关性分析",
            filename = paste0(savepath,"/说明.txt"))
    
    return()
    
  }
  
  corrdata <- corrcal(x = x,
                      y = y,
                      ci = ci,
                      xname = xname,
                      yname = yname,
                      corfilter = corfilter,
                      corfiltertype = corfiltertype,
                      pfilter = pfilter,
                      adjust = adjust,
                      transx = transx,
                      transy = transy,
                      method = cormethod,
                      use = use)
  
  if(is.null(corrdata)){
    
    savetxt(data = "数据为空,不进行相关性分析",
            filename = paste0(savepath,"/说明.txt"))

    return(corrdata)
    
  }
  
  plotdata <- mapmoudle(corr = corrdata$rawcor$r,
                        p.mat = corrdata$rawcor$p.adj,
                        mapname = mapname,
                        savepath = savepath,
                        ...)
  
  corrdata$plotdata <- plotdata
  
  if (savecorr) {
    # savexlsx1(data = corrdata$filterdata$linkdata, 
    #           filename = paste0(savepath,"/",mapname,".xlsx"), 
    #           sheet = "cordata")
    
    savetxt(data = corrdata$framedata$linkdata,
            filename = paste0(savepath,"/",mapname,".xls"))
  } 
  
  return(corrdata)
}

