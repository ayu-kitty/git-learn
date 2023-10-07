
#' @export
dataID_old <- function(code.negID,code.posID,...){
  print("--------------------")
  print("--------------------")
  print("定性数据处理")
  print("---负离子数据获取---")
  negID <- IDfetch_old(name=code.negID)
  print("---正离子数据获取---")
  posID <- IDfetch_old(name=code.posID)
  print("---正/负离子合并---")
  negposID <- negposIDmerge_old(datanegID=negID,dataposID=posID)
  print("---定性数据筛选---")
  negposID <- datascreen_old(data=negposID,...)
  #negposID <- classified(data=negposID)
  print("定性数据处理完成")
  print("--------------------")
  print("--------------------")
  return(negposID)
}


#ID文件获取
#' @export
IDfetch_old <- function(name){
  tryfetch <- try({
    print("定性数据读取")
    data <- readdata(filename,errorstop = T)
    print("定性数据提取")
    data <- collectlist(data,needlist = c("Compound","Description","Compound ID",
                                          "Adducts","Formula","Score",
                                          "Fragmentation Score","Mass Error (ppm)","Charge"))
    names(data)[1] <- "ID"
    names(data)[2] <- "Metabolites"
    names(data)[9] <- "Ion mode"
    
    # print("定性数据筛选")
    # data <- datascreen(data=data,...)
  },silent = F)
  
  if(class(tryfetch)=="try-error"){
    data <- NULL
    warning("数据获取失败",immediate. = T)
  }
  return(data)
}

#negID和posID合并
#' @export
negposIDmerge_old <- function(datanegID,dataposID){
  ####ID文件数据去重####
  negID<-datanegID
  posID<-dataposID
  if(!is.null(negID)){negID$`Ion mode`<-"neg"}
  if(!is.null(posID)){posID$`Ion mode`<-"pos"}
  ####ID文件数据neg和pos合并，去重####
  negposID<-rbind(negID,posID)#将-neg-ID.csv和-pos-ID.csv数据合并
  #negposID<-simpledatascreen(data=negposID)
  return(negposID)
}


#用于定性文件的筛选
#' @export
datascreen_old <- function(data,
                       peptide = F,
                       meta = NULL,
                       org = NULL,
                       score = 36,
                       fragscore = 0){
  predata <- data
  
  if (!peptide) {
    print("定性数据肽段剔除，无需剔除添加参数peptide=T")
    range <- as.character(c(7641:24040,64961:65360,103480:265032))
    predata <- data[!(data$`Compound ID` %in% range),]
  }
  
  print(paste0("score<",score,"剔除"))
  predata <- predata[!predata$Score<score, ]
  print(paste0("Fragmentation Score<",fragscore,"剔除"))
  predata <- predata[!predata$`Fragmentation Score`<fragscore, ]
  
  # 定性结果处理
  predata <- predata[predata$Metabolites!="", ]
  predata<-tidyr::separate(data = predata, col ="Metabolites", into = c("Metabolites"), sep = " /",extra = "drop")
  
  print("添加自建库level-***")
  predata[,"level"]<- "***"
  predata[grep(pattern = "^LM",x = predata$`Compound ID`),"level"] <- ""
  predata[grep(pattern = "^HMDB",x = predata$`Compound ID`),"level"] <- ""
  predata[grep(pattern = "^[0-9]*$",x = predata$`Compound ID`),"level"] <- ""
  predata$name <- tolower(predata$Metabolites)
  
  if (!is.null(meta)) {
    print("特选代谢物加分筛选")
    meta1 <- tolower(meta)
    predata[predata$name %in% meta1,"Score"] <- predata[predata$name %in% meta1,"Score"]+100
  }
  
  print("添加分类及kegg")
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
  predata[is.na(predata$`Super Class`),"Super Class"] <- "Unclassified"
  predata[is.na(predata$`Class`),"Class"] <- "Unclassified"
  predata[is.na(predata$`Sub Class`),"Sub Class"] <- "Unclassified"
  print("添加kegg库level-**及添加人动值库level-*")
  if(!is.null(org)){
    if(org == "map"){
      predata[(!is.na(predata$kegg))& predata$level == "","level"] <- "**"
    }else{
      METAdata <- getdatabase(databasewd = databasepath(path = "kegg"),database = "keggbackground.db",name = org)
      predata[(predata$kegg %in%  METAdata$KEGG) & predata$level == "","level"] <- "**"
      METAdata <- getdatabase(database = "other.db",name = "HMDB")
      if(org == "hsa"){
        predata[(predata$id %in% METAdata[METAdata$Endogenous == 1, "id"]) & predata$level == "","level"] <- "*"
      }else{
        METAdata2 <- getdatabase(databasewd = databasepath(path = "kegg"),database = "kegg.db",name = "organism")
        org2 <- METAdata2[METAdata2$abbreviation == org,"class2"]
        if(length(org2)==0){
        }else if(org2 == "Animals"){predata[(predata$id %in%  METAdata[METAdata$Animal == 1, "id"]) & predata$level == "","level"] <- "*"
        }else if(org2 == "Plants"){predata[(predata$id %in%  METAdata[METAdata$Plant == 1, "id"]) & predata$level == "","level"] <- "*"}
      }
    }
  }else{predata[(!is.na(predata$kegg))& predata$level == "","level"] <- "**"}
  
  #predata <- data1
  print("去除数据重复")
  predata1<-predata
  
  while(dim(predata)[1]!=0){
    predata1<-predata1[order(predata1$level,decreasing = T),]
    predata1<-predata1[order(predata1$Score,decreasing = T),]
    predata1<-predata1[!duplicated(predata1$ID), ]
    
    #predata1<-predata1[!duplicated(predata1$id), ]
    predata1<-predata1[!duplicated(predata1$name), ]
    #predata1<-predata1[!duplicated(predata1$ID), ]
    
    predata<-predata[!(predata$ID %in% predata1$ID),]
    #predata<-predata[!(predata$id %in% predata1$id),]
    predata<-predata[!(predata$name %in% predata1$name),]
    predata1<-rbind(predata1,predata)
  }
  predata<-predata1
  colnames(predata)[20]<-"kegg"
  predata <- predata[,c(-1)]
  predata <- predata[,c("ID","Ion mode",
                        "Metabolites","Compound ID",
                        "Super Class","Class","Sub Class",
                        "level","kegg",
                        "Score","Fragmentation Score",
                        "Adducts","Formula",
                        "Mass Error (ppm)")]
  
  return(predata)
}
