#!/opt/conda/bin/Rscript

#' @export
get_alldata_file <- function(datafile = "oecloud/rawdata/datafile.txt",
                             infofile = "oecloud/rawdata/infofile.txt",
                             savepath = NULL){
  
  # 数据提取
  data <- readdata(datafile,row.names = 1)
  info <- readdata(infofile,row.names = 1)
  # info <- info[,1,drop = F]
  data <- merge(info,data,by = 0,sort = F)[,-1,drop = F]
  
  if(is.null(savepath)){
    
    return(data)
    
  }else{
    # 数据保存
    savetxt(data = data,
            filename = savepath)
  }
  
}

#' @export
get_alldata_file_cloud <- function(datafile = "oecloud/rawdata/datafile.txt",
                                   infofile = "oecloud/rawdata/infofile.txt",
                                   classfile = "oecloud/rawdata/classfile.yaml",
                                   comparefile = "oecloud/rawdata/compare.yaml",
                                   savepath = NULL){
  
  # 数据提取
  data <- readdata(datafile,row.names = 1)
  info <- readdata(infofile,row.names = 1)
  # info <- info[,1,drop = F]
  data <- merge(info,data,by = 0,sort = F)[,-1,drop = F]
  
  class <- readdata(classfile)
  # samplename <- unique(unlist(class))
  # class2 <- data.frame("Sample" = samplename,row.names = samplename)
  # groupname <- names(class)
  # for ( i in 1:length(groupname)) {
  #   class2[class[[groupname[i]]],groupname[i]] <- groupname[i]
  # }
  class2 <- NULL
  for (i in 1:length(class)) {
    class3 <- data.frame("Sample" = class[[i]],"Group" = names(class)[i])
    class2 <- rbind(class2,class3)
  }
  class <- class2
  class <- class[class[,1] %in% colnames(data),]
  
  compare <- readdata(comparefile)
  compare2 <- names(compare$params$all_compare)
  if(!is.null(compare2)){
    compare2 <- data.frame("比较组" = compare2)
    for ( i in 1:dim(compare2)[1]) {
      compare2[i,"配对"] <- compare$compare[[compare2[i,"比较组"]]]$paired
    }
    compare2[,"比较组"] <- gsub(pattern = "-vs-",replacement = "/",x = compare2[,"比较组"])
  }
  compare <- compare2
  
  if(is.null(savepath)){
    
    return(list(data = data,
                class = class,
                compare = compare))
    
  }else{
    
    wb <- openxlsx::createWorkbook()
    addsheet1(data = data,wb = wb,sheet = "数据矩阵")
    addsheet1(data = class,wb = wb,sheet = "分组")
    addsheet1(data = compare,wb = wb,sheet = "比较")
    savewb(wb = wb,filename = savepath)
    
  }
  
}

#' @export
get_alldata_file_cloud2 <- function(rawdatafile = "oecloud/rawdata/rawdatafile.txt",
                                    datafile = "oecloud/rawdata/datafile.txt",
                                    infofile = "oecloud/rawdata/infofile.txt",
                                    classfile = "oecloud/rawdata/classfile.yaml",
                                    comparefile = "oecloud/rawdata/compare.yaml",
                                    savepath = NULL){
  
  # 数据提取

  
  data <- readdata(datafile,row.names = 1)
  info <- readdata(infofile,row.names = 1)
  # info <- info[,1,drop = F]
  data <- merge(info,data,by = 0,sort = F)[,-1,drop = F]
  
  if(file.exists(rawdatafile)){
    rawdata <- readdata(rawdatafile,row.names = 1)
    rawdata <- merge(info,rawdata,by = 0,sort = F)[,-1,drop = F]
  }else{
    rawdata <- NULL
  }
  
  class <- readdata(classfile)
  # samplename <- unique(unlist(class))
  # class2 <- data.frame("Sample" = samplename,row.names = samplename)
  # groupname <- names(class)
  # for ( i in 1:length(groupname)) {
  #   class2[class[[groupname[i]]],groupname[i]] <- groupname[i]
  # }
  class2 <- NULL
  for (i in 1:length(class)) {
    class3 <- data.frame("Sample" = class[[i]],"Group" = names(class)[i])
    class2 <- rbind(class2,class3)
  }
  class <- class2
  class <- class[class[,1] %in% colnames(data),]
  
  compare <- readdata(comparefile)
  compare2 <- names(compare$params$all_compare)
  if(!is.null(compare2)){
    compare2 <- data.frame("比较组" = compare2)
    for ( i in 1:dim(compare2)[1]) {
      compare2[i,"配对"] <- compare$compare[[compare2[i,"比较组"]]]$paired
    }
    compare2[,"比较组"] <- gsub(pattern = "-vs-",replacement = "/",x = compare2[,"比较组"])
  }
  compare <- compare2
  
  if(is.null(savepath)){
    
    return(list(rawdata = rawdata,
                data = data,
                class = class,
                compare = compare))
    
  }else{
    
    wb <- openxlsx::createWorkbook()
    addsheet1(data = rawdata,wb = wb,sheet = "缺失值数据矩阵")
    addsheet1(data = data,wb = wb,sheet = "数据矩阵")
    addsheet1(data = class,wb = wb,sheet = "分组")
    addsheet1(data = compare,wb = wb,sheet = "比较")
    savewb(wb = wb,filename = savepath)
    
  }
  
}
