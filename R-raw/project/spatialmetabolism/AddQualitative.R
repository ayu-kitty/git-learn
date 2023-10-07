#' AddQualitative
#'
#' 在空代定性的基础上以mz添加定性结果
#'
#' @param mzdata mz数据
#' @param qualitativefrom 定性数据库位置
#' @param database 使用数据库类型
#' @param ppm mz偏差
#' @param meta 关注代谢物
#'
#' @return 定性矩阵
#'
#' @export
AddQualitative <- function(mzdata,
                           mode,
                           ppm = 10,
                           meta = NA,
                           addmzinfo = F) {
  
  mz <- mzdata
  
  # 获取代谢物
  acumeta <- getmysqldata(dbname = "meta",table = "hmdbmzforsptial",wherename = "mode",wheredata = mode)
  
  # hmdbdata <- getmysqldata(table = "hmdbmzforsptial2",
  #                          wherename = "Adduct",
  #                          wheredata = c("M-","M-2H","M+2H","M+NH4"),
  #                          wherenotin = T)
  # delcpdid <- getLCMSUntargetedDelCpdsforhmdbid()
  # hmdbdata <- hmdbdata[!(hmdbdata$`Compound ID` %in% delcpdid),]
  # hmdbdata <- hmdbdata[!grepl(pattern = "K",x = hmdbdata$Formula),]
  # hmdbdata <- hmdbdata[!grepl(pattern = "Na",x = hmdbdata$Formula),]
  # hmdbdata <- hmdbdata[!grepl(pattern = "Li",x = hmdbdata$Formula),]
  # hmdbdata <- hmdbdata[!grepl(pattern = "Ca",x = hmdbdata$Formula),]
  # hmdbdata <- hmdbdata[!grepl(pattern = "Cu",x = hmdbdata$Formula),]
  # hmdbdata <- hmdbdata[!grepl(pattern = "Mg",x = hmdbdata$Formula),]
  # hmdbdata <- hmdbdata[!grepl(pattern = "Co",x = hmdbdata$Formula),]
  # hmdbdata <- hmdbdata[!grepl(pattern = "Cr",x = hmdbdata$Formula),]
  # hmdbdata <- hmdbdata[!grepl(pattern = "Fe",x = hmdbdata$Formula),]
  # hmdbdata <- hmdbdata[grepl(pattern = "^C",x = hmdbdata$Formula),]
  # acumeta <- hmdbdata[hmdbdata$mode == mode,]
  
  acumeta <- acumeta[order(acumeta$`Compound ID`), ]
  
  newname <- getmysqldata(dbname = "meta",table = "hmdbnewname")
  acumeta <- merge(acumeta,newname,by.x= "Compound ID",by.y = "hmdbid")
  acumeta$Metabolites <- acumeta$Name
  acumeta <- acumeta[,colnames(acumeta) != "Name"]
  acumeta <- acumeta[,colnames(acumeta) != "mode"]
  acumeta <- acumeta[!duplicated(acumeta[,c("Metabolites","Adduct")]), ]
  
  acumeta <- acumeta[!is.na(acumeta$mz),]
  acumeta[,"Adduct"] <- gsub(pattern = ".*M",replacement = "",x = acumeta[,"Adduct"])
  acumeta[,"moc"] <- NA
  acumeta[,"spec"] <- NA
  acumeta[,"spat"] <- NA
  acumeta[,"msm"] <- NA
  acumeta[,"project"] <- NA
  
  matchmz2 <- data.frame()
  # 匹配
  for (i in 1:length(mz)) {
    minmz <- mz[i] - ppm / 1000000 * mz[i]
    maxmz <- mz[i] + ppm / 1000000 * mz[i]
    matchmz <- acumeta[((acumeta$mz <= maxmz) & (acumeta$mz >= minmz)), ]
    if (dim(matchmz)[1] > 0) {
      
      matchmz <- matchmz[order(matchmz$mz),]
      if (any(matchmz$`Compound ID` %in% meta)) {
        print(paste0("mz:", mz[i, "mz"], "以特选代谢物定性"))
        matchmz <- matchmz[matchmz$id %in% meta, ]
      }
      
      if(addmzinfo){
        matchmz[,"realmz"] <- mz[i]
      }
      
      matchmz2 <- rbind(matchmz2,matchmz)
    }
  }
  
  return(matchmz2)
}



#' manQualitative
#'
#' 搭建人工数据库
#'
#' @param hmdbid hmdbid编号
#'
#' @export
manQualitative <- function(hmdbid) {
  
  # 获取代谢物
  acumeta <- getmysqldata(dbname = "meta",table = "hmdbmzforsptial",wherename = "`Compound ID`",wheredata = hmdbid)
  acumeta <- acumeta[order(acumeta$`Compound ID`), ]
  newname <- getmysqldata(dbname = "meta",table = "hmdbnewname")
  acumeta <- merge(acumeta,newname,by.x= "Compound ID",by.y = "hmdbid")
  acumeta$Metabolites <- acumeta$Name
  acumeta <- acumeta[,colnames(acumeta) != "Name"]
  acumeta <- acumeta[!duplicated(acumeta[,c("Metabolites","Adduct")]), ]
  
  acumeta <- acumeta[!is.na(acumeta$mz),]
  acumeta[,"Adduct"] <- gsub(pattern = ".*M",replacement = "",x = acumeta[,"Adduct"])
  acumeta[,"moc"] <- NA
  acumeta[,"spec"] <- NA
  acumeta[,"spat"] <- NA
  acumeta[,"msm"] <- NA
  acumeta[,"project"] <- NA
  acumeta[,"level"] <- "***"
  
  savexlsx1(data = acumeta[acumeta$mode == "neg",colnames(acumeta) != "mode"],name = "人工数据库.xlsx",sheet = "neg")
  savexlsx1(data = acumeta[acumeta$mode == "pos",colnames(acumeta) != "mode"],name = "人工数据库.xlsx",sheet = "pos")
}
