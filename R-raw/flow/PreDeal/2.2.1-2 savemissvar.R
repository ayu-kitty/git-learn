#!/opt/conda/bin/Rscript


#' 保存原始数据
#' 
#' @param filename 读取文件名
#' @param savefilename 保存文件名
#' @param method 缺失值处理方法
#'
#' @export
missvalueprocess <- function(filename = "data_raw.csv",
                             savefilename = "data_missimpute.csv",
                             method = "halfmin"){
  
  print("缺失值填充")
  data <- readdata(filename = filename,row.names = 1)
  data[is.na(data)] <- 0 
  dataline <- as.numeric(matrix(as.matrix(data[,-1]),ncol = 1))
  dataline <- dataline[!is.na(dataline)]
  dataline <- dataline[dataline!=0]
  if(any(dataline < 0) & method == "halfmin"){
    warning("由于矩阵中有小数，最小值一半填充替换为最小值填充",immediate. = T)
    method <- "valuemin"
  }
  
  if(method == "halfmin"){
    data[,-1][data[,-1] == 0] <- min(dataline)/2
  }else if(method == "valuemin"){
    data[,-1][data[,-1] == 0] <- min(dataline)
  }
  
  write.csv(x = data,file = savefilename,row.names = T)
  return(savefilename)
}


#' missvaluedel
#' 
#' 保存原始数据
#' 
#' @param filename 读取文件名
#' @param savefilename 保存文件名
#' @param method 缺失值处理方法
#'
#' @export
missvaluedel <- function(filename = "data_raw.csv",
                         savefilename = "data_missdel.csv",
                         missvarremovepercent = 0.5,
                         missvarremovebygroup = T){
  
  print("缺失值删除")
  data <- readdata(filename = filename,row.names = 1)
  data[is.na(data)] <- 0 

  if(missvarremovebygroup){
    
    data <- data[,c(T,apply(data[,-1], 2, mulmissvalueratio,group = data[,1]) <= missvarremovepercent)]
      
  }else{
    
    data <- data[,c(T,apply(data[,-1], 2, missvalueratio) <= missvarremovepercent)]
  
  }
  
  write.csv(x = data,file = savefilename,row.names = T)
  return(savefilename)
}



#' @export
missvalueratio <- function(x){
  sum(x==0)/length(x)
}

#' @export
mulmissvalueratio <- function(x,group){
  
  singlegroup <- unique(group)
  rationum <- NULL
  
  for ( i in singlegroup) {
    
    data <- x[group == i]
    rationum2 <- missvalueratio(data)
    rationum <- c(rationum,rationum2)
    
  }
  
  return(min(rationum))
  
}


