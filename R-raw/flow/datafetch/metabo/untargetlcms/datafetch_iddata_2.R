#!/opt/conda/bin/Rscript

#' 获取定性数据
#'
#' @param code.negID 负离子定性路径
#' @param code.posID 正离子定性路径
#' @param ... 见`datascreen()`
#'
#' @export
dataID_2 <- function(code.negID, 
                     code.posID, 
                     ...) {
  print("--------------------")
  print("--------------------")
  print("定性数据处理")
  
  print("---负离子数据获取---")
  negID <- IDfetch(name = code.negID)
  
  print("---正离子数据获取---")
  posID <- IDfetch(name = code.posID)
  
  print("---正/负离子合并---")
  negposID <- negposIDmerge(datanegID = negID, 
                            dataposID = posID)
  
  print("---定性数据筛选---")
  negposID <- datascreen_2(data = negposID, ...)
  
  print("定性数据处理完成")
  print("--------------------")
  print("--------------------")
  return(negposID)
}

#' 用于定性文件的筛选
#'
#' @param data 定性数据
#' @param peptide 逻辑，是否保留肽段
#' @param meta 关注代谢物
#' @param score score筛选标准
#' @param fragscore 二级筛选标准
#'
#' @export
datascreen_2 <- function(data,
                         peptide = F,
                         meta = NULL,
                         score = 36,
                         fragscore = 0,
                         Drug = T,
                         Tcm = F,
                         ...) {
  predata <- data
  
  if (!peptide) {
    print("~定性数据肽段剔除，无需剔除添加参数peptide=T")
    range <- as.character(c(7641:24040, 64961:65360, 103480:265032))
    predata <- predata[!(predata$`Compound ID` %in% range), ]
  }
  
  print(paste0("score<", score, "剔除"))
  predata <- predata[!predata$Score < score, ]
  print(paste0("Fragmentation Score<", fragscore, "剔除"))
  predata <- predata[!predata$`Fragmentation Score` < fragscore, ]
  
  # 定性结果处理
  predata <- predata[predata$Metabolites != "", ]
  predata <- getmetainfo(data = predata,idlist = "Compound ID")
  predata <- predata[!is.na(predata$cid), ]
  
  predata[,"level"] <- "Level 4"
  predata[(predata[,"Fragmentation Score"] > 45),"level"] <- "Level 3"
  predata[,"rtscore"] <- predata$Score - (predata$`Fragmentation Score`/5)-40
  predata[(predata[,"rtscore"] > 0),"level"] <- "Level 2"
  predata[(predata[,"rtscore"] > 0) & (predata[,"Fragmentation Score"] > 45),"level"] <- "Level 1"
  
  if (!is.null(meta)) {
    print("关注代谢物加****")
    meta1 <- tolower(meta)
    predata[tolower(predata$Metabolites) %in% meta1, "level"] <- "Level 0"
  }
  
  if (Drug) {
    print("~正在进行Drugblank去除,不进行这项,请使用Drug = F")
    delcpds  <- getLCMSUntargetedDelCpds()[2]
    predata <- predata[!(predata$InChIKey %in% delcpds), ]
  }else{
    print("~不进行Drugblank去除,如果要去除,请使用Drug = T")
  }
  
  print("去除数据重复")
  predata <- predata[order(predata$Score, decreasing = T), ]
  predata <- predata[order(predata$level), ]
  predata <- predata[!duplicated(predata$ID), ]
  predata <- predata[!duplicated(predata$cid), ]
  predata <- predata[!duplicated(tolower(predata$`Name`)), ]
  predata <- predata[!is.na(predata$`Name`), ]
  
  if(Tcm){
    predata <- predata[, c(
      "ID", "Ion mode",
      "Name","cid", "HMDB","METLIN","Lipidmaps","KEGG","ChEBI","PubChem","CAS","smiles","InChIKey",
      "Super Class", "Class", "Sub Class","中文名","中文大类","中文子类",
      "level","Score", "Fragmentation Score",
      "Adducts", "Formula",
      "Mass Error (ppm)"#,"厂家","货号","纯度"
    )]
  }else{
    predata <- predata[, c(
      "ID", "Ion mode",
      "Name","cid", "HMDB","METLIN","Lipidmaps","KEGG","ChEBI","PubChem","CAS","smiles","InChIKey",
      "Super Class", "Class", "Sub Class",
      "level","Score", "Fragmentation Score",
      "Adducts", "Formula",
      "Mass Error (ppm)"#,"厂家","货号","纯度"
    )]
    
  }

  colnames(predata)[colnames(predata) == "Name"] <- "Metabolites"
  
  selectshortname <- function(x){
    name <- unlist(strsplit(x,split = ";"))
    b1 <- stringr::str_count(string = name,pattern = "\\(")
    b2 <- stringr::str_count(string = name,pattern = "\\)")
    if(all(b1 == b2)){
      name <- name[which.min(nchar(name))]
    }else{
      name <- x
    }
    return(name)
  }
  
  predata$Metabolites <- apply(predata[,"Metabolites",drop = F], 1,selectshortname)
  predata <-  predata[nchar(predata$Metabolites) < 100,]
    
  return(predata)
}

