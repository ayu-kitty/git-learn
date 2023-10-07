#!/opt/conda/bin/Rscript

# setwd("/data/hstore1/database/test/2023-04-21差异分析测试")

#' 整理数据
#' 
#' @param rawfile 原始文件名称
#' @param datafile 数据保存路径
#' @param infofile 信息保存路径
#' @param classfile 分组保存路径
#'
#' @export
organizedata <- function(rawfile = "数据矩阵.xlsx",
                         rawdatafile = c(rawfile,"数据矩阵"),
                         rawclassfile = c(rawfile,"分组"),
                         rawcomparefile = c(rawfile,"比较"),
                         datafile = "oecloud/rawdata/datafile.txt",
                         infofile = "oecloud/rawdata/infofile.txt",
                         inputdatafile = paste0(dirname(datafile),"/inputfile.txt"),
                         classfile = "classfile.yaml",
                         classtypefile = "classtype.xlsx",
                         comparefile = "oecloud/rawdata/compare.yaml",
                         compare = NULL,
                         p_method = "t.test",
                         mul_p_method = "oneway.test",
                         p_adjust_method = "BH",
                         paired = F,
                         keylist = NULL,
                         overwrite = F,
                         ...){
  if(is.character(rawdatafile)){if(rawdatafile[1] == ""){rawdatafile <- c(rawfile,"数据矩阵")}}
  if(is.character(rawclassfile)){if(rawclassfile[1] == ""){rawclassfile <- c(rawfile,"分组")}}
  if(is.character(rawcomparefile)){if(rawcomparefile[1] == ""){rawcomparefile <- c(rawfile,"比较")}}
  
  paired <- as.logical(paired)
  if(is.na(paired)){paired <- F}
  
  # 数据读取
  rawdata <- readdata(filename = rawdatafile)
  
  if(!is.data.frame(rawdata) & is.list(rawdata)){
    
    organizeobj(obj = rawdata,
                datafile = datafile,
                infofile = infofile ,
                classfile = classfile,
                classtypefile = classtypefile,
                comparefile = comparefile,
                p_method = p_method,
                mul_p_method = mul_p_method,
                p_adjust_method = p_adjust_method,
                paired = paired,
                overwrite = overwrite,
                ...)
    
    return()
  }
  
  savetxt(data = rawdata,filename = inputdatafile,overwrite = overwrite)
  
  classtry <- try({
    class <- readdata(filename = rawclassfile,row.names = 1)
  })
  
  if("try-error" %in% class(classtry)){
    class <- data.frame()
  }
  
  if(dim(class)[1] == 0){
    class <- readdata(filename = rawclassfile)
    samplename <- unique(class[,1])
    class2 <- data.frame("sample" = samplename,row.names = samplename)
    class2 <- class2[,0,drop = F]
    groupname <- unique(class[,2])
    for ( i in 1:length(groupname)) {
      class2[class[class[,2] == groupname[i],1],groupname[i]] <- groupname[i]
    }
    class <- class2
  }
  
  # if(any(is.na(class[,1]))){
  #   
  #   if(dim(class)[2] > 1){
  #     for ( i in 2:dim(class)[2]) {
  #       
  #     }
  #   }
  #   
  # }
  
  # [x] 由于第一列有重复列名报错，使用自定义列名
  # row.names(rawdata) <- rawdata[,1]
  if(!is.null(keylist)){
    if(max(table(rawdata[,keylist])) == 1){
      row.names(rawdata) <- rawdata[,keylist]
    }else{
      warning("指定行名有重复使用自定义行名")
      row.names(rawdata) <- paste0("_",1:dim(rawdata)[1])
    }
  }else{
    row.names(rawdata) <- paste0("_",1:dim(rawdata)[1])
  }
  
  rawdata <- rawdata[,!(colnames(rawdata) %in% row.names(class)[apply(class, 1, function(x){ any( x== "删除",na.rm = T)})]),drop = F]
  class <- class[!apply(class, 1, function(x){ any( x== "删除",na.rm = T)}),,drop = F]
  class <- class[,!apply(class,2, function(x){all(is.na(x))}),drop = F]
  
  # 数据提取
  if(!all(rownames(class) %in% colnames(rawdata))){
    noclass <- class[!(rownames(class) %in% colnames(rawdata)),,drop = F]
    noclass <-  paste0(rownames(noclass),collapse = ";")
    stop(paste0(noclass,"样本不在数据矩阵中"))
  }
  
  data <- rawdata[,colnames(rawdata) %in% rownames(class),drop = F]
  
  # 数据信息提取
  infodata <- rawdata[,!(colnames(rawdata) %in% rownames(class)),drop = F]
  
  # 分组信息提取
  newclass <- list()
  for ( i in 1:dim(class)[2]) {
    newclass2 <- class[,i,drop = F]
    newclass2 <- newclass2[!is.na(newclass2[,1]),,drop = F]
    newclass2 <- newclass2[newclass2[,1]!="",,drop = F]
    if(dim(newclass2)[1] == 0){
      next
    }
    uniqueclass <- unique(newclass2[,1])
    if("QC" %in% uniqueclass){
      uniqueclass <- uniqueclass[uniqueclass!="QC"]
      uniqueclass <- c("QC",uniqueclass)
    }
    for( j in 1:length(uniqueclass)){
      if(uniqueclass[j] %in% names(newclass)){
        stop(paste0(uniqueclass[j],"为重复分组名"))
      }
      if(!grepl(pattern = "^[A-Za-z0-9._+-]*$",x = uniqueclass[j])){
        stop(paste0(uniqueclass[j],"分组有特殊字符,请勿使用;/,，等特殊字符"))
      }
      newclass3 <- newclass2[newclass2[,1] == uniqueclass[j],,drop = F]
      newclass4 <- list(row.names(newclass3))
      names(newclass4) <- uniqueclass[j]
      newclass <- c(newclass,newclass4)
    }
  }
  
  main_group <- NULL
  if(length(newclass) > 2){
    main_group <- names(newclass)[1]
    group_sample <- newclass[[1]]
    for ( i in 2:length(newclass)) {
      if(!any(newclass[[i]] %in% group_sample)){
        main_group <- c(main_group,names(newclass)[i])
        group_sample <- c(group_sample,newclass[[i]])
      }
    }
  }else{
    main_group <- names(newclass)
  }
  
  if("QC" %in% main_group){
    main_group <- main_group[main_group!="QC"]
    main_group <- c("QC",main_group)
  }
  main_group2 <- main_group[main_group!="QC"]
  
  config <- list(params = list("all_group" = names(newclass),
                               "main_group" = main_group,
                               "main_group-qc" = main_group2),
                 compare = list("All" = list(name = paste(main_group,collapse = "-vs-"),
                                             rawname = paste(main_group,collapse = "/"),
                                             group = main_group,
                                             groupnum = length(main_group),
                                             samplenum = sum(lengths(newclass[main_group])),
                                             minsamplenum = min(lengths(newclass[main_group])),
                                             maxsamplenum = max(lengths(newclass[main_group])),
                                             paired = F,
                                             "p_method"= ifelse(length(main_group) > 2,mul_p_method,p_method),
                                             "p_adjust_method" = p_adjust_method),
                                "Allsample" = list(name = paste(main_group2,collapse = "-vs-"),
                                                   rawname = paste(main_group2,collapse = "/"),
                                                   group = main_group2,
                                                   groupnum = length(main_group2),
                                                   samplenum = sum(lengths(newclass[main_group2])),
                                                   minsamplenum = min(lengths(newclass[main_group2])),
                                                   maxsamplenum = max(lengths(newclass[main_group])),
                                                   paired = F,
                                                   "p_method"= ifelse(length(main_group2) > 2,mul_p_method,p_method),
                                                   "p_adjust_method" = p_adjust_method)))
  
  classtype <- data.frame(Group = names(newclass))
  classtype[,"fill"] <- stylefun_group(classfile = classtype,
                                       styletype = "fill",
                                       ...)
  classtype[,"colour"] <- stylefun_group(classfile = classtype,
                                         styletype = "colour",
                                         ...)
  if(dim(classtype)[1] > 1){
    if(classtype$Group[1] == "QC"){
      classtype[2:dim(classtype)[1],"fill"] <- classtype[1:(dim(classtype)[1]-1),"fill"]
      classtype[1,"fill"] <- "#A6D854"
      if(classtype[1,"colour"] != classtype[2,"colour"]){
        classtype[2:dim(classtype)[1],"colour"] <- classtype[1:(dim(classtype)[1]-1),"colour"]
        classtype[1,"colour"] <- "#A6D854"
      }
    }
  }

  classtype[,"shape"] <- stylefun_group(classfile = classtype,
                                        styletype = "shape",
                                        ...)
  
  saveorganizedata(data = data,
                   infodata = infodata,
                   classdata = newclass,
                   classtype = classtype,
                   datafile = datafile,
                   infofile = infofile,
                   classfile = classfile,
                   classtypefile = classtypefile,
                   overwrite = overwrite)
  
  comparetry <- try({
    if(!is.null(compare)){
      compare <- data.frame(group = compare)
    }else{
      compare <- readdata(filename = rawcomparefile)
      compare <- compare[!is.na(compare$比较组),,drop = F]
    }
  })
  
  if("try-error" %in% class(comparetry)){
    saveyaml(data = config,filename = comparefile,overwrite = overwrite)
  }else if(dim(compare)[1] == 0){
    saveyaml(data = config,filename = comparefile,overwrite = overwrite)
  }else{
    
    colnames(compare)[colnames(compare) == "比较组"] <- "group"
    
    if(!("配对" %in% colnames(compare))){
      if("paired" %in% colnames(compare)){
        
      }else{
        compare[,"paired"] <- paired
      }
    }else{
      colnames(compare)[colnames(compare) == "配对"] <- "paired"
    }
    
    compare[,"paired"] <- as.logical(compare[,"paired"])
    compare[is.na(compare[,"paired"]),"paired"] <- paired
    
    config[["params"]][["all_compare"]] <- list()
    
    for ( i in 1:dim(compare)[1]) {
      main_group <- compare[i,"group"]
      main_group <- unlist(strsplit(main_group,split = "/"))
      
      if(!all(main_group %in% names(newclass))){
        stop(paste0(paste0(main_group[!(main_group %in% names(newclass))],collapse = ";")," 不在分组中"))
      }
      
      compare2 <- list(list(name = paste(main_group,collapse = "-vs-"),
                            rawname = paste(main_group,collapse = "/"),
                            group = main_group,
                            groupnum = length(main_group),
                            samplenum = sum(lengths(newclass[main_group])),
                            minsamplenum = min(lengths(newclass[main_group])),
                            maxsamplenum = max(lengths(newclass[main_group])),
                            paired = compare[i,"paired"],
                            "p_method"= ifelse(length(main_group) > 2,mul_p_method,p_method),
                            "p_adjust_method" = p_adjust_method))
      
      if(length(main_group) > 2){compare2[[1]][["paired"]] <- F}
      if(length(main_group) == 2){
        if(lengths(newclass[main_group[1]]) != lengths(newclass[main_group[2]])){
          compare2[[1]][["paired"]] <- F
        }
      }
      
      names(compare2) <- paste(main_group,collapse = "-vs-")
      config[["params"]][["all_compare"]] <- c(config[["params"]][["all_compare"]],compare2)
    }
    
    config[["compare"]] <- c(config[["compare"]],config[["params"]][["all_compare"]])
    
    saveyaml(data = config,filename = comparefile,overwrite = overwrite)
    
  }
  
}


#' 对于meta包中的obj进行整理数据
#' 
#' @param obj 数据
#' @param datafile 数据保存路径
#' @param infofile 信息保存路径
#' @param classfile 分组保存路径
#'
#' @export
organizeobj <- function(obj,
                        datafile = "oecloud/rawdata/datafile.txt",
                        infofile = "oecloud/rawdata/infofile.txt",
                        classfile = "classfile.yaml",
                        classtypefile = "classtype.xlsx",
                        comparefile = "oecloud/rawdata/compare.yaml",
                        compare = NULL,
                        p_method = "t.test",
                        mul_p_method = "oneway.test",
                        p_adjust_method = "BH",
                        paired = F,
                        overwrite = F,
                        ...){
  
  obj <- readdata(filename = obj)
  
  if(any(is.character(obj))){
    obj <- get(obj[1])
  }
  
  newclass <- strsplit(obj[["info"]][["class"]][["samples"]],split = ",")
  names(newclass) <- obj[["info"]][["class"]][["group"]]
  
  for( j in 1:length(newclass)){
    if(!grepl(pattern = "^[A-Za-z0-9._+-]*$",x =names(newclass)[j])){
      stop(paste0(names(newclass)[j],"分组有特殊字符,请勿使用;/,，等特殊字符"))
    }
  }
  
  main_group <- unique(obj[["info"]][["sample"]][["class"]][[1]])
  if("QC" %in% main_group){
    main_group <- main_group[main_group!="QC"]
    main_group <- c("QC",main_group)
  }
  main_group2 <- main_group[main_group!="QC"]
  
  config <- list(params = list("all_group" = names(newclass),
                               "main_group" = main_group,
                               "main_group-qc" = main_group2),
                 compare = list("All" = list(name = paste(main_group,collapse = "-vs-"),
                                             rawname = paste(main_group,collapse = "/"),
                                             group = main_group,
                                             groupnum = length(main_group),
                                             samplenum = sum(lengths(newclass[main_group])),
                                             minsamplenum = min(lengths(newclass[main_group])),
                                             maxsamplenum = max(lengths(newclass[main_group])),
                                             paired = F,
                                             "p_method"= ifelse(length(main_group) > 2,mul_p_method,p_method),
                                             "p_adjust_method" = p_adjust_method),
                                "Allsample" = list(name = paste(main_group2,collapse = "-vs-"),
                                                   rawname = paste(main_group2,collapse = "/"),
                                                   group = main_group2,
                                                   groupnum = length(main_group2),
                                                   samplenum = sum(lengths(newclass[main_group2])),
                                                   minsamplenum = min(lengths(newclass[main_group2])),
                                                   maxsamplenum = max(lengths(newclass[main_group])),
                                                   paired = F,
                                                   "p_method"= ifelse(length(main_group2) > 2,mul_p_method,p_method),
                                                   "p_adjust_method" = p_adjust_method)))
  
  classtype <- data.frame(Group = names(newclass))
  classtype[,"fill"] <- stylefun_group(classfile = classtype,
                                       styletype = "fill",
                                       ...)
  classtype[,"colour"] <- stylefun_group(classfile = classtype,
                                         styletype = "colour",
                                         ...)
  if(dim(classtype)[1] > 1){
    if(classtype$Group[1] == "QC"){
      classtype[2:dim(classtype)[1],"fill"] <- classtype[1:(dim(classtype)[1]-1),"fill"]
      classtype[1,"fill"] <- "#A6D854"
        if(classtype[1,"colour"] != classtype[2,"colour"]){
          classtype[2:dim(classtype)[1],"colour"] <- classtype[1:(dim(classtype)[1]-1),"colour"]
          classtype[1,"colour"] <- "#A6D854"
        }
    }
  }
  
  classtype[,"shape"] <- stylefun_group(classfile = classtype,
                                        styletype = "shape",
                                        ...)
  
  if(!is.null(obj[["data"]][["zeroprocess"]][["data"]])){
    prerawdatafile <- paste0(dirname(datafile),"/zerodatafile.txt")
    if(!file.exists(prerawdatafile)){
      savetxt(data = obj[["data"]][["zeroprocess"]][["data"]],filename = prerawdatafile,row.names = T)
    }
  }
  
  saveorganizedata(data = obj[["data"]][["data"]][["data"]],
                   infodata = obj[["data"]][["predata"]][["information"]],
                   classdata = newclass,
                   classtype = classtype,
                   datafile = datafile,
                   infofile = infofile,
                   classfile = classfile,
                   classtypefile = classtypefile,
                   overwrite = overwrite)
  
  comparetry <- try({
    if(!is.null(compare)){
      compare <- data.frame(group = compare)
    }else{
      compare <- data.frame(group = unlist(obj$info$dif$compare),
                            paired = unlist(obj$info$dif$paired))
    }
  })
  
  if("try-error" %in% class(comparetry)){
    saveyaml(data = config,filename = comparefile,overwrite = overwrite)
  }else if(dim(compare)[1] == 0){
    saveyaml(data = config,filename = comparefile,overwrite = overwrite)
  }else{
    
    colnames(compare)[colnames(compare) == "比较组"] <- "group"
    
    if(!("配对" %in% colnames(compare))){
      if("paired" %in% colnames(compare)){
        
      }else{
        compare[,"paired"] <- paired
      }
    }else{
      colnames(compare)[colnames(compare) == "配对"] <- "paired"
    }
    
    compare[,"paired"] <- as.logical(compare[,"paired"])
    compare[is.na(compare[,"paired"]),"paired"] <- paired
    
    config[["params"]][["all_compare"]] <- list()
    
    for ( i in 1:dim(compare)[1]) {
      main_group <- compare[i,"group"]
      main_group <- unlist(strsplit(main_group,split = "/"))
      
      if(!all(main_group %in% names(newclass))){
        stop(paste0(paste0(main_group[!(main_group %in% names(newclass))],collapse = ";")," 不在分组中"))
      }
      
      compare2 <- list(list(name = paste(main_group,collapse = "-vs-"),
                            rawname = paste(main_group,collapse = "/"),
                            group = main_group,
                            groupnum = length(main_group),
                            samplenum = sum(lengths(newclass[main_group])),
                            minsamplenum = min(lengths(newclass[main_group])),
                            maxsamplenum = max(lengths(newclass[main_group])),
                            paired = compare[i,"paired"],
                            "p_method"= ifelse(length(main_group) > 2,mul_p_method,p_method),
                            "p_adjust_method" = p_adjust_method))
      
      if(length(main_group) > 2){compare2[[1]][["paired"]] <- F}
      if(length(main_group) == 2){
        if(lengths(newclass[main_group[1]]) != lengths(newclass[main_group[2]])){
          compare2[[1]][["paired"]] <- F
        }
      }
      
      names(compare2) <- paste(main_group,collapse = "-vs-")
      config[["params"]][["all_compare"]] <- c(config[["params"]][["all_compare"]],compare2)
    }
    
    config[["compare"]] <- c(config[["compare"]],config[["params"]][["all_compare"]])
    
    saveyaml(data = config,filename = comparefile,overwrite = overwrite)
    
  }
}

#' @export
saveorganizedata <- function(data,
                             infodata,
                             classdata,
                             classtype,
                             datafile = "oecloud/rawdata/datafile.txt",
                             infofile = "oecloud/rawdata/infofile.txt",
                             classfile = "classfile.yaml",
                             classtypefile = "classtype.xlsx",
                             overwrite = F,
                             ...){
  
  if(!all(apply(data,2,is.numeric))){
    stop(paste0(paste(colnames(data)[!apply(data,2,is.numeric)],collapse = ";"),"列有非数值文本"))
  }
  
  # 数据保存
  savetxt(data = data,filename = datafile,row.names = T,overwrite = overwrite)
  savetxt(data = infodata,filename = infofile,row.names = T,overwrite = overwrite)
  saveyaml(data = classdata,filename = classfile,overwrite = overwrite)
  savexlsx1(data = classtype,filename = classtypefile,overwrite = overwrite)
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-rf","--rawfile",default = "数据矩阵.xlsx", help = "数据矩阵")
  parser$add_argument("-rdf","--rawdatafile",default = NULL, help = "表达数据")
  parser$add_argument("-rcf","--rawclassfile",default = NULL, help = "分组数据")
  parser$add_argument("-rcpf","--rawcomparefile",default = NULL, help = "比较组数据")
  parser$add_argument("-fp","--fillpalette",default = NULL, help = "填充颜色pattle")
  parser$add_argument("-fop","--filloutpalette",default = NULL, help = "外部填充颜色pattle",nargs = "+")
  parser$add_argument("-cp","--colourpalette",default = NULL, help = "边框颜色pattle")
  parser$add_argument("-cop","--colouroutpalette",default = NULL, help = "外部边框颜色pattle",nargs = "+")
  
  parser$add_argument("-df","--datafile",default = "oecloud/rawdata/datafile.txt", help = "数据保存路径")
  parser$add_argument("-if","--infofile",default = "oecloud/rawdata/infofile.txt", help = "信息保存路径")
  parser$add_argument("-cf","--classfile",default = "classfile.yaml", help = "分组保存路径")
  parser$add_argument("-ctf","--classtypefile",default = "classtype.xlsx", help = "分组绘图保存路径")
  parser$add_argument("-cpf","--comparefile",default = "oecloud/rawdata/compare.yaml", help = "比较组保存路径")
  
  args <- parser$parse_args()
  
  args <- args[!unlist(lapply(args,is.null))]
  
  flowargs <- do.call(what = organizedata,args = args)
  
}
