#!/opt/conda/bin/Rscript

#' 获取定性数据
#'
#' @param code.negID 负离子定性路径
#' @param code.posID 正离子定性路径
#' @param ... 见`datascreen()`
#'
#' @export
dataID <- function(code.negID, 
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
  negposID <- datascreen(data = negposID, ...)
  
  print("定性数据处理完成")
  print("--------------------")
  print("--------------------")
  return(negposID)
}

#' ID文件获取
#'
#' @param name LCMS的ID文件名
#'
#' @export
IDfetch <- function(name) {
  tryfetch <- try(
    {
      print("定性数据读取")
      data <- readdata(filename = name,filetype = "csv")
      
      if(colnames(data)[1] == ""){
        data <- readdata(filename = name,filetype = "csv",skip = 1)
      }
      
      print("定性数据提取")
      needlist <- c("Compound", "Description", "Compound ID",
                    "Adducts", "Formula", "Score",
                    "Fragmentation Score", "Mass Error (ppm)", "Charge",
                    "Retention time (min)")
      data <- data[,needlist,drop = F]
      names(data)[1] <- "ID"
      names(data)[2] <- "Metabolites"
      names(data)[9] <- "Ion mode"
      
    },
    silent = F
  )
  
  if (class(tryfetch) == "try-error") {
    data <- NULL
    warning("数据获取失败", immediate. = T)
  }
  return(data)
}

#' negID和posID合并
#'
#' @param datanegID 负离子ID数据
#' @param dataposID 正离子ID数据
#'
#' @export
negposIDmerge <- function(datanegID, 
                          dataposID) {
  #### ID文件数据去重####
  negID <- datanegID
  posID <- dataposID
  
  if (!is.null(negID)) {
    negID$`Ion mode` <- "neg"
  }
  if (!is.null(posID)) {
    posID$`Ion mode` <- "pos"
  }
  
  # 将-neg-ID.csv和-pos-ID.csv数据合并
  negposID <- rbind(negID, posID) 
  
  return(negposID)
}

#' 获取RT数据文件
#'
#' @param filename 数据文件名
#'
#' @export
getRTdata <- function(path = lmbio::databasepath(path = "database/qualitative/LCMS")){
  
  filename <- list.files(path = path,pattern = "^RT",full.names = T)
  filename <- filename[order(filename,decreasing = T)]
  
  RTdata <- readdata(filename[1])
  
  RTdata <- RTdata[,c("inchiKey","Retention time (min)","source")]
  colnames(RTdata) <- c("InChIKey","RealRT","RTsource")
  return(RTdata)
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
datascreen <- function(data,
                       RTdata = getRTdata(),
                       peptide = F,
                       meta = NULL,
                       score = 36,
                       fragscore = 0,
                       retentiontime = 0.2,
                       adjustRT = T,
                       Drug = T,
                       c17realrentionime = 11.17,
                       clrealrentionime = 3.8,
                       c17diffrentionime = 0.2,
                       cldiffrentionime = 0.15) {
  predata <- data
  
  if (!peptide) {
    print("定性数据肽段剔除，无需剔除添加参数peptide=T")
    range <- as.character(c(7641:24040, 64961:65360, 103480:265032))
    predata <- predata[!(predata$`Compound ID` %in% range), ]
  }
  
  print(paste0("score<", score, "剔除"))
  predata <- predata[!predata$Score < score, ]
  print(paste0("Fragmentation Score<", fragscore, "剔除"))
  predata <- predata[!predata$`Fragmentation Score` < fragscore, ]
  
  # 定性结果处理
  predata <- predata[predata$Metabolites != "", ]
  predata <- tidyr::separate(data = predata, col = "Metabolites", into = c("Metabolites"), sep = " /", extra = "drop")
  predata$name <- tolower(predata$Metabolites)
  
  iddata_mysql <- getmysqldata(dbname = "meta",
                               table = "lmcidtometaid",
                               wherename = "id",
                               wheredata = predata$`Compound ID`)
  iddata_mysql <- iddata_mysql[iddata_mysql$DatabaseFrom %in% c("HMDB","METLIN","Lipidmaps"),]
  iddata_mysql <- iddata_mysql[iddata_mysql$Facticity == 1,]
  iddata_mysql <- iddata_mysql[!duplicated(iddata_mysql[,c("id","DatabaseFrom")]),]
  iddata_mysql <- iddata_mysql[,c("lmcid","id")]
  
  predata <- merge(x = predata,y = iddata_mysql,by.x = "Compound ID",by.y="id",all.x = T)
  
  predata2 <- predata[is.na(predata$lmcid),]
  predata <- predata[!is.na(predata$lmcid),]
  iddata_mysql <- getmysqldata(dbname = "meta",
                               table = "lmcinfo",
                               tablelist = c("lmcid,InChIKey"),
                               wherename = "InChIKey",
                               wheredata = predata$`Compound ID`)
  predata2 <- predata2[,colnames(predata2) != "lmcid"]
  predata2 <- merge(x = predata2,y = iddata_mysql,by.x = "Compound ID",by.y="InChIKey",all.x = T)
  predata <- rbind(predata,predata2)
  
  predata2 <- predata[is.na(predata$lmcid),]
  predata <- predata[!is.na(predata$lmcid),]
  iddata_mysql <- getmysqldata(dbname = "meta",
                               table = "lmcname")
  iddata_mysql <- iddata_mysql[iddata_mysql$Name %in% predata2$name,]
  iddata_mysql <- iddata_mysql[!duplicated(iddata_mysql[,c("Name")]),]
  predata2 <- predata2[,colnames(predata2) != "lmcid"]
  predata2 <- merge(x = predata2,y = iddata_mysql,by.x = "name",by.y="Name",all.x = T)
  predata <- rbind(predata,predata2)
  
  lmcid <- predata$lmcid
  lmcid <- lmcid[!is.na(lmcid)]
  lmcid <- lmcid[!duplicated(lmcid)]
  iddata_mysql <- getmysqldata(dbname = "meta",
                               table = "lmcinfo",
                               tablelist = c("lmcid,InChIKey,Name,SuperClass as 'Super Class',Class,SubClass as 'Sub Class'"),
                               wherename = "lmcid",
                               wheredata = lmcid)
  iddata_mysql[is.na(iddata_mysql$`Super Class`),"Super Class"] <- "Unclassified"
  iddata_mysql[is.na(iddata_mysql$`Class`),"Class"] <- "Unclassified"
  iddata_mysql[is.na(iddata_mysql$`Sub Class`),"Sub Class"] <- "Unclassified"
  predata <- merge(x = predata,y = iddata_mysql,by.x = "lmcid",by.y="lmcid",all.x = T)
  
  iddata_mysql <- getmysqldata(dbname = "meta",
                               table = "lmcidtometaid",
                               wherename = "lmcid",
                               wheredata = lmcid)
  iddata_mysql <- iddata_mysql[iddata_mysql$DatabaseFrom %in% c("HMDB","METLIN","Lipidmaps","PubChem","KEGG","CAS"),]
  iddata_mysql <- iddata_mysql[iddata_mysql$Facticity == 1,]
  iddata_mysql <- iddata_mysql[!duplicated(iddata_mysql[,c("lmcid","DatabaseFrom")]),]
  iddata_mysql <- iddata_mysql[,c("lmcid","id","DatabaseFrom")]
  iddata_mysql <- reshape2::dcast(iddata_mysql,lmcid~DatabaseFrom,value.var = "id")
  predata <- merge(x = predata,y = iddata_mysql,by.x = "lmcid",by.y="lmcid",all.x = T)
  
  if(adjustRT){
    RTdata2 <- RTdata[is.na(RTdata$RTsource),]
    predata2 <- merge(x = predata,y = RTdata2,by.x = "InChIKey",by.y="InChIKey")
    predata2 <- predata2[order(predata2$`Fragmentation Score`,decreasing = T),]
    predata2 <- predata2[order(predata2$`Score`,decreasing = T),]
    predata2 <- predata2[!duplicated(predata2[,c("InChIKey","ID")]),]
    predata2[,"diftime"] <- abs(predata2[,"Retention time (min)"]-predata2[,"RealRT"])
    predata2[,"diftimem"] <- predata2[,"diftime"]/predata2[,"RealRT"]
    predata2[,"diftime2"] <- predata2[,"Retention time (min)"]-predata2[,"RealRT"]
    predata2 <- predata2[predata2$diftimem < 0.15,]
    predata2 <- predata2[predata2$diftime < 1,]
    predata2 <- predata2[predata2$`Fragmentation Score` > 50,]
    predata2 <- predata2[predata2$Score > 50,]
    predata2 <- predata2[predata2$`Retention time (min)`> 0.5,]
    delmeta <- names(table(predata2$InChIKey)[table(predata2$InChIKey) > 1])
    predata2 <- predata2[!(predata2$InChIKey %in% delmeta),]
    
    predata3 <- predata[!is.na(predata$InChIKey),]
    predata3 <- predata3[predata3$InChIKey == "SRRQPVVYXBTRQK-XMMPIXPASA-N",]
    predata3 <- predata3[order(predata3$`Fragmentation Score`,decreasing = T),]
    predata3 <- predata3[order(predata3$`Score`,decreasing = T),]
    predata3 <- predata3[!duplicated(predata3[,c("InChIKey","ID")]),]
    predata3 <- predata3[predata3$Adducts %in% c("M+FA-H","M+H, M+K, M+Na","M+H"),]
    predata3 <- predata3[predata3$`Retention time (min)` > 10.5 & predata3$`Retention time (min)` < 12,]
    if(dim(predata3)[1] > 0){
      c17rentionime <- predata3$`Retention time (min)`[1]
    }else{
      c17rentionime <- NULL
    }
    
    predata4 <- predata[!is.na(predata$InChIKey),]
    predata4 <- predata4[predata4$InChIKey == "XCKUCSJNPYTMAE-UHFFFAOYSA-N",]
    predata4 <- predata4[order(predata4$`Fragmentation Score`,decreasing = T),]
    predata4 <- predata4[order(predata4$`Score`,decreasing = T),]
    predata4 <- predata4[!duplicated(predata4[,c("InChIKey","ID")]),]
    predata4 <- predata4[predata4$Adducts %in% c("M-H","M+H, M+K, M+Na","M+H"),]
    predata4 <- predata4[predata4$`Retention time (min)` > 3 & predata4$`Retention time (min)` < 4.5,]
    if(dim(predata4)[1] > 0){
      clrentionime <- predata4$`Retention time (min)`[1]
    }else{
      clrentionime <- NULL
    }
    
    adjust <- T
    if(!is.null(c17rentionime) & !is.null(clrentionime)){
      print(paste0("C17保留时间：",c17rentionime))
      print(paste0("2-氯苯丙氨酸保留时间：",clrentionime))
      
      if((c17rentionime-c17realrentionime > c17diffrentionime) & (clrentionime-clrealrentionime > cldiffrentionime)){
        predata2 <- predata2[predata2$diftime2 > 0,]
      }else if((c17rentionime-c17realrentionime < -c17diffrentionime) & (clrentionime-clrealrentionime < -cldiffrentionime)){
        predata2 <- predata2[predata2$diftime2 < 0,]
      }else if((abs(c17rentionime-c17realrentionime) < c17diffrentionime) & ( abs(clrentionime-clrealrentionime) < cldiffrentionime)){
        adjust <- F
        print("保留时间偏差不大，不进行校正")
      }else{
        stop("内标的保留时间往不同方向偏移")
      }
    }else if(!is.null(c17rentionime)){
      print(paste0("C17保留时间：",c17rentionime))
      
      if((c17rentionime-c17realrentionime > c17diffrentionime)){
        predata2 <- predata2[predata2$diftime2 > 0,]
      }else if((c17rentionime-c17realrentionime < -c17diffrentionime)){
        predata2 <- predata2[predata2$diftime2 < 0,]
      }else if((abs(c17rentionime-c17realrentionime) < c17diffrentionime)){
        adjust <- F
        print("保留时间偏差不大，不进行校正")
      }
    }else if(!is.null(clrentionime)){
      print(paste0("2-氯苯丙氨酸保留时间：",clrentionime))
      
      if((clrentionime-clrealrentionime > cldiffrentionime)){
        predata2 <- predata2[predata2$diftime2 > 0,]
      }else if((clrentionime-clrealrentionime < -cldiffrentionime)){
        predata2 <- predata2[predata2$diftime2 < 0,]
      }else if((abs(clrentionime-clrealrentionime) < cldiffrentionime)){
        adjust <- F
        print("保留时间偏差不大，不进行校正")
      }
    }else{
      stop("未找到内标")
    }
    
    if(adjust){
      if(dim(predata2)[2] > 10){
        lmdata <- lmcal(x = predata2$RealRT,
                        y = predata2$`Retention time (min)`,
                        clrealrentionime = clrealrentionime,
                        c17realrentionime = c17realrentionime,
                        clrentionime = clrentionime,
                        c17rentionime = c17rentionime,
                        cldiffrentionime = cldiffrentionime,
                        c17diffrentionime = c17diffrentionime)
        summarydata <- summary(lmdata)
        rtmap(lmdata,summarydata)
        RTdata$RealRT <- predict(lmdata,data.frame(x=RTdata$RealRT))
      }else{
        stop("构建模型数据过少，不进行校正")
      }
    }
  }
  
  predata <- merge(x = predata,y = RTdata,by.x = "InChIKey",by.y="InChIKey",all.x = T)
  predata[,"diftime"] <- abs(predata[,"Retention time (min)"]-predata[,"RealRT"])
  predata[is.na(predata[,"diftime"]),"diftime"] <- 10
  predata[,"newScore"] <- predata$Score - predata$`Fragmentation Score`/5
  predata[,"timeScore"] <- round((retentiontime-predata$diftime)/retentiontime*20)
  predata[predata$timeScore < 0,"timeScore"] <- 0
  predata[predata$newScore < 40,"Score"] <- predata[predata$newScore < 40,"Score"]+predata[predata$newScore < 40,"timeScore"]
  
  predata[,"level"] <- "*"
  predata[(predata[,"diftime"] < retentiontime),"level"] <- "**"
  predata[(predata[,"diftime"] < retentiontime) & (predata[,"Fragmentation Score"] > 30),"level"] <- "***"
  predata <- predata[!is.na(predata$lmcid),]
  
  if (!is.null(meta)) {
    print("关注代谢物加****")
    meta1 <- tolower(meta)
    predata[predata$name %in% meta1, "level"] <- "****"
  }
  
  if (Drug) {
  print("正在进行Drugblank去除,不进行这项,请使用Drug = F")
  delcpds  <- lmbio::getLCMSUntargetedDelCpds()[2]
  predata <- predata[!(predata$InChIKey %in% delcpds), ]
  }else{
    print("不进行Drugblank去除,如果要去除,请使用Drug = T")
  }



  print("去除数据重复")
  predata <- predata[order(predata$Score, decreasing = T), ]
  predata <- predata[order(predata$level, decreasing = T), ]
  predata <- predata[!duplicated(predata$ID), ]
  predata <- predata[!duplicated(predata$lmcid), ]
  predata <- predata[!duplicated(predata$`Name`), ]
  
  # print(colnames(predata))
  predata <- predata[, c(
    "ID", "Ion mode",
    "Name", "HMDB","METLIN","Lipidmaps","KEGG","PubChem","CAS","InChIKey",
    "Super Class", "Class", "Sub Class",
    "level","Score", "Fragmentation Score",
    "Adducts", "Formula",
    "Mass Error (ppm)"
  )]
  
  colnames(predata) <- c(
    "ID", "Ion mode",
    "Metabolites", "HMDB","METLIN","Lipidmaps","kegg","PubChem","CAS","InChIKey",
    "Super Class", "Class", "Sub Class",
    "level","Score", "Fragmentation Score",
    "Adducts", "Formula",
    "Mass Error (ppm)"
  )
  
  return(predata)
}

