#!/opt/conda/bin/Rscript

#' 保存预处理的原始数据
#' 
#' @param data 数据
#' @param class 样本分组
#' @param trans 逻辑值，是否转置
#' @param filename 保存文件名
#'
#' @export
saverawdata <- function(data,
                        class,
                        trans = T,
                        filename = "data_raw.csv"){
  
  if(trans){
    data <- data.frame(t(data),
                       check.names = F,
                       stringsAsFactors = F)
  }
  
  if(is.data.frame(class)){
    class <- class[!is.na(class[,1]),,drop=F]
    
    if(!all(rownames(data) %in% rownames(class))){
      if(all(!(rownames(data) %in% rownames(class)))){
        stop("矩阵数据在分组数据中均未找到对应样本")
      }
      
      noexistsample <- rownames(data)[!(rownames(data) %in% rownames(class))]
      warning(paste("以下样本：",
                    paste(noexistsample,collapse = ","),
                    "在分组数据中不存在"),
              immediate. = T)
    }
    
    data <- data[rownames(data) %in% rownames(class),]
    if(dim(data)[1] < 2){
      stop("矩阵数据中样本过少")
    }
    
    class <- class[rownames(data),1,drop = F]
    colnames(class)[1] <- "Group"
    
  }else{
    stop("class参数传输错误")
  }
  
  rawdata <- cbind(class,data)
  
  if(length(unique(rawdata[,1]))==1){
    warning("仅有一个分组，自动变成两组,请谨慎进行标准化",immediate. = T)
    rawdata2 <- rawdata
    rawdata2[,1] <- paste0(rawdata2[,1],"_new_")
    rawdata <- rbind(rawdata,rawdata2)
  }
  
  for (group in unique(rawdata[,1])) {
    rawdata2 <- rawdata[rawdata[,1] == group,,drop = F]
    while(dim(rawdata2)[1] < 3){
      warning(paste0(group,"组样本少于三个，自动拓展本组数据,请谨慎进行标准化"),immediate. = T)
      row.names(rawdata2) <- paste0(row.names(rawdata2),"_new_")
      rawdata <- rbind(rawdata,rawdata2)
      rawdata2 <- rawdata[rawdata[,1] == group,,drop = F]
    }
  }
  
  write.csv(x = rawdata,file = filename,row.names = T)
  
  args <- list(rawdata = rawdata,
               filename = filename,
               data = data,
               class = class)
  
  return(args)
}
