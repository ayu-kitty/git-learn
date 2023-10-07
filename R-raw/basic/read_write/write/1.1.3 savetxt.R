#!/opt/conda/bin/Rscript

#' 保存txt文件
#'
#' @param data dataframe数据
#' @param filename 文件名
#' @param sep 数据分割方式
#' @param row.names 是否保存行名
#'
#' @export
savetxt <- function(data,
                    filename,
                    sep = "\t",
                    row.names = F,
                    fileEncoding = "utf-8",
                    append = T,
                    overwrite = T,
                    ...) {
  
  if(is.null(data)) {return()}
  if(is.data.frame(data)){
    if(dim(data)[1] == 0) {return()}
  }else{
    print(paste0("在", filename, "中保存"))
    
    dirfilename <- createdir(filename = dirname(filename))
    # print(data)
    
    write(x = data,file = filename,append = append,...)
    
    return()
  }
  
  if(row.names != F){
    data2 <- data.frame(FeatureID = rownames(data),row.names = rownames(data))
    if(!is.logical(row.names)){
      names(data2) <- row.names
    }
    data <- cbind(data2,data)
    row.names <- FALSE
  }
  
  print(paste0("在", filename, "中保存"))
  
  dirfilename <- createdir(filename = dirname(filename))
  
  for ( i in 1:dim(data)[2]) {
    
    data[grepl(pattern = "\"",x = data[,i]),i] <- gsub(pattern = "\"",replacement = "",x =data[grepl(pattern = "\"",x = data[,i]),i],)
    
  }
  
  if(!file.exists(filename) | overwrite){
    write.table(x = data,
                file = filename,
                sep = sep,
                row.names = row.names,
                fileEncoding = fileEncoding,
                ...)
    Sys.chmod(paths = filename,mode = "0777",use_umask = F)
  }
  
  # print("保存完毕")
}
