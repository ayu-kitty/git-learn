#!/opt/conda/bin/Rscript

#' 计算相关性
#' 
#' @param x x数据
#' @param y y数据
#' @param ci 是否计算置信区间
#' @param xname x数据名称
#' @param yname y数据名称
#' @param corfilter 相关性筛选标准
#' @param corfiltertype 相关性筛选方法,包含+-,+,-
#' @param pfilter 显著性筛选标准
#' @param adjust 显著性校正算法,包括BH,holm,hochberg,hommel,bonferroni,BY,fdr,none
#' @param transx 是否转置x数据
#' @param transy 是否转置y数据
#' @param method 相关性计算方法,包括pearson,spearman,kendall
#' @param use 相关性计算方法,包括pairwise,complete
#' @param ... 见psych::corr.test
#'
#' @export 
corrcal <- function(x,
                    y = NULL,
                    ci = FALSE,
                    xname = "Featurex",
                    yname = "Featurey",
                    corfilter = 0.95,
                    corfiltertype = "+-",
                    pfilter = 0.05,
                    adjust = "BH",
                    transx = T,
                    transy = T,
                    method = "pearson",
                    use = "pairwise",
                    ...){
  print("相关性计算中")
  
  suppressMessages(library("stringr"))
  
  if(is.null(x)){
    warning("x数据为空,不进行相关性运算",immediate. = T)
    return(NULL)
  }
  
  if(transx){x <- t(x)}
  if(!is.null(y)){if(transy){y <- t(y)}}
  
  if(dim(x)[2] == 0){
    warning("x数据为空,不进行相关性运算",immediate. = T)
    return(NULL)
  }else if(dim(x)[2] == 0){
    warning("x数据为空,不进行相关性运算",immediate. = T)
    return(NULL)
  }
  
  
  for(k in 1:ncol(x)){
    colnames(x)[k] <- unlist(strsplit(split = ";\n",x = colnames(x)[k]))[1]
    colnames(x)[k] <- unlist(strsplit(split = "; ",x = colnames(x)[k]))[1]
    if(str_length(colnames(x)[k])>70){
      newname <- paste0(substring(colnames(x)[k],1,35),"...")
      i <- 2
      while (newname %in% colnames(x)) {
        newname <- paste0(substring(colnames(x)[k],1,35),"...-",i)
        i <- i+1
      }
      colnames(x)[k] <- newname
    }
  }
  
  if(is.null(y)){
    ct <- psych::corr.test(x = x,y = x,ci = ci,
                           adjust = adjust,
                           method = method,
                           use = use,
                           ...)
  }else{
    
    for(k in 1:ncol(y)){
      colnames(y)[k] <- unlist(strsplit(split = ";\n",x = colnames(y)[k]))[1]
      colnames(y)[k] <- unlist(strsplit(split = "; ",x = colnames(y)[k]))[1]
      if(str_length(colnames(y)[k])>70){
        newname <- paste0(substring(colnames(y)[k],1,35),"...")
        i <- 2
        while (newname %in% colnames(y)) {
          newname <- paste0(substring(colnames(y)[k],1,35),"...-",i)
          i <- i+1
        }
        colnames(y)[k] <- newname
      }
    }
    
    ct <- psych::corr.test(x = x,y = y,ci = ci,
                           adjust = adjust,
                           method = method,
                           use = use,
                           ...)
  }
  
  cordata <- list(rawcor = ct)
  
  rdata <- as.data.frame(ct$r)
  pdata <- as.data.frame(ct$p.adj)
  rdata[, xname] <- row.names(rdata)
  pdata[, xname] <- row.names(pdata)
  rdata <- reshape2::melt(rdata, id.vars = xname, variable.name = yname, value.name = "cor", factorsAsStrings = F)
  pdata <- reshape2::melt(pdata, id.vars = xname, variable.name = yname, value.name = "p", factorsAsStrings = F)
  combinedata <- merge(rdata, pdata,by=c(xname,yname),)
  combinedata[,yname] <- as.character(combinedata[,yname])
  
  if(is.null(y)){
    combinedata <- combinedata[combinedata[,xname] != combinedata[,yname],]
    combinedata <- combinedata[!duplicated(apply(X = combinedata[,c(xname,yname)], MARGIN = 1,FUN = orderpaste)),]
    nodedata <- data.frame("name" = colnames(x),
                           "type" = xname)
  }else{
    nodedatax <- data.frame("name" = colnames(x),
                            "type" = xname)
    nodedatay <- data.frame("name" = colnames(y),
                            "type" = yname)
    nodedata <- rbind(nodedatax,nodedatay)
  }
  combinedata[, "class"] <- "positive correlation"
  combinedata[combinedata[, "cor"] < 0, "class"] <- "negative correlation"
  
  cordata$filterdata$corfilter <- corfilter
  cordata$filterdata$pfilter <- pfilter
  
  cordata$framedata$linkdata <- combinedata
  cordata$framedata$nodedata <- nodedata
  
  # 筛选
  if(corfiltertype == "+"){
    filterlinkdata <- combinedata[combinedata$cor > corfilter,]
  }else if(corfiltertype == "-"){
    filterlinkdata <- combinedata[combinedata$cor < -corfilter,]
  }else{
    filterlinkdata <- combinedata[combinedata$cor > corfilter | combinedata$cor < -corfilter,]
  }
  filterlinkdata <- filterlinkdata[filterlinkdata$p < pfilter,]
  
  
  if(dim(filterlinkdata)[1] > 0){
    filternodedata <- nodedata
    nodename <- c(filterlinkdata[,xname],filterlinkdata[,yname])
    nodename <- data.frame(table(nodename))
    filternodedata <- merge(filternodedata,nodename,by.x ="name",by.y = "nodename")
  }else{
    filternodedata <- nodedata
    filternodedata[,"Degree"] <- 0
    filternodedata <- filternodedata[0,]
    warning("相关性筛选后无符合条件的结果",immediate. = T)
  }
  
  cordata$filterdata$linkdata <- filterlinkdata
  cordata$filterdata$nodedata <- filternodedata
  
  # print("相关性计算完成")
  
  return(cordata)
}

#' @export 
orderpaste <- function(x){
  x <- x[order(x)]
  chr <- paste(x,collapse = "-")
  return(chr)
}
