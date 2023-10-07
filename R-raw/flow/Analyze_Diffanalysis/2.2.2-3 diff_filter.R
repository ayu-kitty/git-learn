#!/opt/conda/bin/Rscript

#' @export
all_diff_filter <- function(rdspath = "oecloud/diff_ana/",
                            rdsname = list.files(path = rdspath,pattern = "\\.rds$",full.names = T),
                            fcfilter = c(1.2,1.5,2),
                            noprint = F,
                            ...){
  
  if (length(fcfilter) > 1) {
    
    result <- NULL
    
    for ( i in 1:length(fcfilter)) {
      result2 <- all_diff_filter(rdsname = rdsname,
                                 fcfilter = fcfilter[i],
                                 noprint = T,
                                 savediffrds = NULL,
                                 ...)
      colnames(result2)[2] <-  fcfilter[i]
      if(is.null(result)){
        result <- result2
      }else{
        result <- merge(result,result2,by = "Compare")
      }
    }
    
    print(result)
    
    score <- apply(result[,-1,drop = F], 2, calfitdiffscore)
    fitfc <- as.numeric(names(which.max(score)))
    
    print(paste0("#最佳的FC筛选标准为",fitfc))
    
    savetxt(data = fitfc,filename = paste0(rdspath,"/fitfc.txt"))
    
    result <- all_diff_filter(rdsname = rdsname,
                              fcfilter = fitfc,
                              noprint = noprint,
                              ...)
    
  }else{
    result <- lapply(rdsname,diff_filter,
                     fcfilter = fcfilter,
                     noprint = noprint,
                     ...)
    result <- data.frame(Compare = gsub(pattern = "\\..*",
                                        replacement = "",
                                        x = basename(rdsname)),
                         filternum = unlist(result))
  }
  
  return(result)
}

#' 对筛选后数量进行打分
#' 
#' @export
calfitdiffscore <- function(x,
                            fitnum = 200,
                            adjust = 150){
  if(any(x <= 10)){
    x[x <= 10] <- -100000000
  }
  
  score <- abs(x-fitnum)+1
  # score <- (log2(adjust)-log2(score))/log2(adjust)
  score <- (adjust-score)/adjust
  score <- mean(score)
  
  return(score)
  
}

#' 差异筛选
#' 
#' @export
diff_filter <- function(rdsname,
                        savepath = "oecloud/diff_filter",
                        savediffrds = paste0(savepath,"/diff_filter-",gsub(pattern = "diff_ana-",replacement = "",basename(rdsname))),
                        ...){
  
  data <- readdata(rdsname)
  
  data <- diff_filter_cal(data,
                          savediffrds = savediffrds,
                          ...)
  
  return(data)
}

#' @export
diff_filter_cal <- function(data,
                            savediffrds = "diffdata.rds",
                            vipfilter = 0,
                            pfilter = 0.05,
                            adjpfilter = 0,
                            fcfilter = 0,
                            fcfiltertype ="+-",
                            errorptofcfilter = 2,
                            needlist = c("Metabolites","Accession"),
                            noprint = F){
  
  diffdata <- data$result$diffdata
  infodata <- data$args$info
  
  if("VIP" %in% colnames(diffdata)){
    if(vipfilter > 0){
      diffdata <- diffdata[!is.na(diffdata$VIP),]
      diffdata <- diffdata[diffdata$VIP > vipfilter,]
    }
  }else if(vipfilter > 0){
    vipfilter <- 0
    if(errorptofcfilter > 0 & fcfilter <= 0){
      fcfilter <- errorptofcfilter
    }
  }else{
    vipfilter <- 0
  }
  
  if("p-value" %in% colnames(diffdata)){
    if(pfilter > 0){
      diffdata <- diffdata[!is.na(diffdata$`p-value`),]
      diffdata <- diffdata[diffdata$`p-value` < pfilter,]
    }
    if(adjpfilter > 0){
      pfilter <- 0
      diffdata <- diffdata[!is.na(diffdata$`q-value`),]
      diffdata <- diffdata[diffdata$`q-value` < adjpfilter,]
    }else{
      adjpfilter <- 0
    }
  }else if(pfilter > 0 | adjpfilter > 0){
    pfilter <- 0
    adjpfilter <- 0
    if(errorptofcfilter > 0 & fcfilter <= 0){
      fcfilter <- errorptofcfilter
    }
  }else{
    pfilter <- 0
    adjpfilter <- 0
  }
  
  logfcname <- "log2FoldChange"
  
  if(fcfilter > 0 & logfcname %in% colnames(diffdata)){
    diffdata <- diffdata[!is.na(diffdata[,logfcname]),]
    if(fcfiltertype == "+-"){
      diffdata <- diffdata[diffdata[,logfcname] > log2(fcfilter) | diffdata[,logfcname] < -log2(fcfilter),]
    }else if(fcfiltertype == "+"){
      diffdata <- diffdata[diffdata[,logfcname] > log2(fcfilter),]
    }else if(fcfiltertype == "-"){
      diffdata <- diffdata[diffdata[,logfcname] < -log2(fcfilter),]
    }
  }else{
    fcfilter <- 0
  }
  
  if(!is.null(needlist)){
    
    for ( i in 1:length(needlist)) {
      if(needlist[i] %in% colnames(infodata)){
        infodata <- infodata[!is.na(infodata[,needlist[i]]),]
      }
    }
    
    diffdata <- diffdata[row.names(diffdata) %in% row.names(infodata),]
  }
  
  data$result$filterdata <- diffdata
  data$result$filterargs <- list(vipfilter = vipfilter,
                                 pfilter = pfilter,
                                 adjpfilter = adjpfilter,
                                 fcfilter = fcfilter,
                                 fcfiltertype = fcfiltertype)
  
  class(data) <- c(class(data),"diff_filter")
  
  if(!is.null(savediffrds)){
    saverds(data = data,filename = savediffrds)
  }
  
  if(!noprint){
    print(paste0("#",data$args$compare,"组筛选后有",dim(diffdata)[1],"个"))
  }
  
  return(dim(diffdata)[1])
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-rp","--rdspath",default = "oecloud/diff_ana/", help = "rds数据路径")
  parser$add_argument("-r","--rdsname",
                      default = NULL, 
                      help = "中间过程数据位置",
                      nargs = "+")
  parser$add_argument("-ff","--fcfilter",default = c(1.2,1.5,2), help = "FC筛选标准",type = "double",nargs = "+")
  parser$add_argument("-sp","--savepath",default = "oecloud/diff_filter", help = "保存路径")
  parser$add_argument("-fft","--fcfiltertype",default = "+-", help = "FC筛选标准")
  parser$add_argument("-eff","--errorptofcfilter",default = 2, help = "FC筛选标准",type = "double")
  parser$add_argument("-vf","--vipfilter",default = 0, help = "vip筛选标准",type = "double")
  parser$add_argument("-pf","--pfilter",default = 0.05, help = "p-value筛选标准",type = "double")
  parser$add_argument("-apf","--adjpfilter",default = 0, help = "adjp-value筛选标准",type = "double")
  
  args <- parser$parse_args()
  
  if(is.null(args$rdsname)){
    args$rdsname <- NULL
  }
  
  flowargs <- do.call(what = all_diff_filter,args = args)
  
}
