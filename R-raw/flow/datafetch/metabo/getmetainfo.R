#!/opt/conda/bin/Rscript

#' @export
getmetainfo <- function(data,
                        databasefrom = c("inchikey","metlin","hmdb","lipidmaps","smiles","cas","kegg","英文名"),
                        idlist = "Compound ID",
                        vague = F,
                        distinct = F,
                        needlist = c("Name","cid",
                                     "CAS","HMDB","METLIN","Lipidmaps","KEGG","PubChem","ChEBI",
                                     "smiles","InChIKey","Mass","Formula",
                                     "Super Class","Class","Sub Class","中文名","中文大类","中文子类",
                                     "厂家","货号","纯度","dbfrom","来源"),
                        cid = F){
  predata <- data[,idlist,drop = F]
  colnames(predata) <- "Compound ID"
  predata <- predata[!is.na(predata$`Compound ID`),,drop = F]
  predata <- predata[predata$`Compound ID`!="",,drop = F]
  predata <- predata[!duplicated(predata$`Compound ID`),,drop = F]
  
  if(cid){
    iddata_mysql <- predata
    iddata_mysql[,"cid"] <- predata$`Compound ID`
  }else{
    iddata_mysql <- NULL
    
    if("inchikey" %in% databasefrom){
      iddata_mysql_inchikey <- getmysqldata(dbname = "cosa",
                                            table = "compound_structure",
                                            tablelist = "cid,inchikey as 'Compound ID'",
                                            wherename = "inchikey",
                                            wheredata = predata$`Compound ID`)
      iddata_mysql <- rbind(iddata_mysql,iddata_mysql_inchikey) 
    }
    
    if("metlin" %in% databasefrom){
      iddata_mysql_metlin <- getmysqldata(dbname = "cosa",
                                          table = "cid_metlinid",
                                          tablelist = "cid,metlinid as 'Compound ID'",
                                          wherename = "metlinid",
                                          wheredata = predata$`Compound ID`)
      iddata_mysql <- rbind(iddata_mysql,iddata_mysql_metlin)
    }
    
    if("hmdb" %in% databasefrom){
      iddata_mysql_hmdb <- getmysqldata(dbname = "cosa",
                                        table = "cid_hmdbid",
                                        tablelist = "cid,hmdbid as 'Compound ID'",
                                        wherename = "hmdbid",
                                        wheredata = predata$`Compound ID`)
      iddata_mysql <- rbind(iddata_mysql,iddata_mysql_hmdb)
    }
    
    if("lipidmaps" %in% databasefrom){
      iddata_mysql_lipid <- getmysqldata(dbname = "cosa",
                                         table = "cid_lipidmapsid",
                                         tablelist = "cid,lipidmapsid as 'Compound ID'",
                                         wherename = "lipidmapsid",
                                         wheredata = predata$`Compound ID`)
      iddata_mysql <- rbind(iddata_mysql,iddata_mysql_lipid)
    }
    
    if("smiles" %in% databasefrom){
      iddata_mysql_smiles <- getmysqldata(dbname = "cosa",
                                          table = "compound_structure",
                                          tablelist = "cid,smiles as 'Compound ID'",
                                          wherename = "smiles",
                                          wheredata = predata$`Compound ID`)
      iddata_mysql <- rbind(iddata_mysql,iddata_mysql_smiles) 
    }
    
    if("kegg" %in% databasefrom){
      iddata_mysql_kegg <- getmysqldata(dbname = "cosa",
                                        table = "cid_keggid",
                                        tablelist = "cid,keggid as 'Compound ID'",
                                        wherename = "keggid",
                                        wheredata = predata$`Compound ID`)
      iddata_mysql <- rbind(iddata_mysql,iddata_mysql_kegg)
    }
    
    if("cas" %in% databasefrom){
      iddata_mysql_syn <- getmysqldata(dbname = "cosa",
                                       table = "compound_identifier",
                                       tablelist = "cid,casNumber as 'Compound ID'",
                                       wherename = "casNumber",
                                       wheredata = predata$`Compound ID`)
      iddata_mysql <- rbind(iddata_mysql,iddata_mysql_syn)
    }
    
    if("英文名" %in% databasefrom){
      if(vague){
        distinct <- T
        iddata_mysql_syn <- NULL
        for ( i in 1:length(predata$`Compound ID`)) {
          iddata_mysql_syn2 <- getmysqldata_vague(dbname = "cosa",
                                                  table = "synonym",
                                                  tablelist = "cid,name as 'Compound ID'",
                                                  wherename = "name",
                                                  wheredata = tolower(predata$`Compound ID`[i]))
          if(dim(iddata_mysql_syn2)[1] == 0){
            iddata_mysql_syn <- rbind(iddata_mysql_syn,iddata_mysql_syn2)
            next
          }
          iddata_mysql_syn2[,'Compound ID'] <- predata$`Compound ID`[i]
          iddata_mysql_syn <- rbind(iddata_mysql_syn,iddata_mysql_syn2)
        }
      }else{
        iddata_mysql_syn <- getmysqldata(dbname = "cosa",
                                         table = "synonym",
                                         tablelist = "cid,name as 'Compound ID'",
                                         wherename = "name",
                                         wheredata = tolower(predata$`Compound ID`))
      }
      iddata_mysql <- rbind(iddata_mysql,iddata_mysql_syn)
    }
    
    if("中文名" %in% databasefrom){
      if(vague){
        distinct <- T
        iddata_mysql_herb <- NULL
        for ( i in 1:length(predata$`Compound ID`)) {
          iddata_mysql_herb2 <- getmysqldata_vague(dbname = "cosa",
                                                   table = "herb_annotation",
                                                   tablelist = "inchikey,cn_name as 'Compound ID'",
                                                   wherename = "cn_name",
                                                   wheredata = predata$`Compound ID`[i])
          if(dim(iddata_mysql_herb2)[1] == 0){
            iddata_mysql_herb <- rbind(iddata_mysql_herb,iddata_mysql_herb2)
            next
          }
          iddata_mysql_herb2[,'Compound ID'] <- predata$`Compound ID`[i]
          iddata_mysql_herb <- rbind(iddata_mysql_herb,iddata_mysql_herb2)
        }
        
      }else{
        iddata_mysql_herb <- getmysqldata(dbname = "cosa",
                                          table = "herb_annotation",
                                          tablelist = "inchikey,cn_name as 'Compound ID'",
                                          wherename = "cn_name",
                                          wheredata = predata$`Compound ID`)
      }
      
      iddata_mysql_inchikey <- getmysqldata(dbname = "cosa",
                                            table = "compound_structure",
                                            tablelist = "cid,inchikey",
                                            wherename = "inchikey",
                                            wheredata = iddata_mysql_herb$inchikey)
      iddata_mysql_herb <- merge(iddata_mysql_inchikey,iddata_mysql_herb,by = "inchikey")[,-1]
      iddata_mysql <- rbind(iddata_mysql,iddata_mysql_herb) 
    }
    
    # base::print(databasefrom)
    if("来源" %in% databasefrom){
      distinct <- T
      iddata_mysql_syn <- getmysqldata_vague(dbname = "cosa",
                                             table = "herb_ac",
                                             tablelist = "cid,Herb_cn_name as 'Compound ID'",
                                             wherename = "Herb_cn_name",
                                             wheredata = predata$`Compound ID`)
      iddata_mysql <- rbind(iddata_mysql,iddata_mysql_syn)
      # base::print(iddata_mysql)
    }
    
    if(!distinct){
      iddata_mysql <- iddata_mysql[!duplicated(iddata_mysql$`Compound ID`),]
    }
    
    if(dim(iddata_mysql)[1] == 0){
      data[,"newlist"] <- tolower(data[,idlist])
      predata <- merge(x = data,y = iddata_mysql,by.x = "newlist",by.y = "Compound ID",all.x = T)[,-1,drop = F]
      return(predata)
    }
    # iddata_mysql[,"cidurl"] <- paste0("http://funmeta.oebiotech.cn/network/",iddata_mysql[,"cid"])
    
  }
  
  iddata_mysql_info <- getmysqldata(dbname = "cosa",
                                    table = "compound_identifier",
                                    tablelist = "cid,name as Name,casNumber as CAS,hmdbid as HMDB,metlinid as METLIN,lipidmapsid as Lipidmaps,keggid as KEGG,pubchemid as PubChem",
                                    wherename = "cid",
                                    wheredata = iddata_mysql$cid)
  iddata_mysql <- merge(iddata_mysql,iddata_mysql_info,by = "cid",all.x = T)
  
  iddata_mysql_inchikey <- getmysqldata(dbname = "cosa",
                                        table = "compound_structure",
                                        tablelist = "cid,inchikey as InChIKey,smiles",
                                        wherename = "cid",
                                        wheredata = iddata_mysql$cid)
  iddata_mysql <- merge(iddata_mysql,iddata_mysql_inchikey,by = "cid",all.x = T)
  
  iddata_mysql_chebi_hmdb <- getmysqldata(dbname = "cosa",
                                          table = "extra_hmdb",
                                          tablelist = "cid,chebiid as ChEBI",
                                          wherename = "cid",
                                          wheredata = iddata_mysql$cid)
  iddata_mysql_chebi_lipid <- getmysqldata(dbname = "cosa",
                                           table = "extra_lipidmaps",
                                           tablelist = "cid,chebiid as ChEBI",
                                           wherename = "cid",
                                           wheredata = iddata_mysql$cid)
  iddata_mysql_chebi <- rbind(iddata_mysql_chebi_hmdb,iddata_mysql_chebi_lipid)
  iddata_mysql_chebi <- iddata_mysql_chebi[!duplicated(iddata_mysql_chebi$cid),]
  iddata_mysql <- merge(iddata_mysql,iddata_mysql_chebi,by = "cid",all.x = T)
  
  iddata_mysql_property <- getmysqldata(dbname = "cosa",
                                        table = "compound_computed_property",
                                        tablelist = "cid,computedExactMass as Mass,computedFormula as Formula",
                                        wherename = "cid",
                                        wheredata = iddata_mysql$cid)
  iddata_mysql_property$Mass <- as.numeric(iddata_mysql_property$Mass)
  iddata_mysql <- merge(iddata_mysql,iddata_mysql_property,by = "cid",all.x = T)
  
  iddata_mysql_classification <- getmysqldata(dbname = "cosa",
                                              table = "compound_classification",
                                              tablelist = "cid,superclass as 'Super Class',class as Class,subclass as 'Sub Class'",
                                              wherename = "cid",
                                              wheredata = iddata_mysql$cid)
  iddata_mysql <- merge(iddata_mysql,iddata_mysql_classification,by = "cid",all.x = T)
  
  iddata_mysql_hreb <- getmysqldata(dbname = "cosa",
                                    table = "herb_annotation",
                                    tablelist = "inchikey as InChIKey,cn_name as '中文名',cn_class1 as '中文大类',cn_class2 as '中文子类',factory as '厂家',item_number as '货号',purity as '纯度',db as dbfrom",
                                    wherename = "inchikey",
                                    wheredata = iddata_mysql$InChIKey)
  iddata_mysql_hreb[is.na(iddata_mysql_hreb$厂家) | iddata_mysql_hreb$厂家== "nan","厂家"] <- NA
  iddata_mysql_hreb[is.na(iddata_mysql_hreb$货号) | iddata_mysql_hreb$货号== "nan","货号"] <- NA
  iddata_mysql_hreb[is.na(iddata_mysql_hreb$纯度) | iddata_mysql_hreb$纯度== "nan","纯度"] <- NA
  iddata_mysql <- merge(iddata_mysql,iddata_mysql_hreb,by = "InChIKey",all.x = T)
  
  iddata_mysql[,"来源"] <- NA
  iddata_mysql_herb <- getmysqldata(dbname = "cosa",
                                    table = "herb_ac",
                                    tablelist = "cid,Herb_cn_name",
                                    wherename = "cid",
                                    wheredata = iddata_mysql$cid)
  for ( i in 1:dim(iddata_mysql)[1]) {
    if(iddata_mysql$cid[i] %in% iddata_mysql_herb$cid){
      if(is.na(iddata_mysql$dbfrom[i])){
        iddata_mysql$dbfrom[i] <- "Herb"
      }else{
        iddata_mysql$dbfrom[i] <- paste0(iddata_mysql$dbfrom[i],";Herb")
      }
      Herb_cn_name <- iddata_mysql_herb[iddata_mysql_herb$cid == iddata_mysql[i,"cid"],"Herb_cn_name"]
      Herb_cn_name <- Herb_cn_name[!duplicated(Herb_cn_name)]
      iddata_mysql[i,"来源"] <- paste(Herb_cn_name,collapse = ";")
    }
  }
  
  iddata_mysql[,"Compound ID"] <- tolower(iddata_mysql[,"Compound ID"])
  iddata_mysql[is.na(iddata_mysql$`Super Class`),"Super Class"] <- "Unclassified"
  iddata_mysql[is.na(iddata_mysql$`Class`),"Class"] <- "Unclassified"
  iddata_mysql[is.na(iddata_mysql$`Sub Class`),"Sub Class"] <- "Unclassified"
  iddata_mysql <- dplyr::select(iddata_mysql,all_of(c("Compound ID",needlist)))
  
  data[,"newlist"] <- tolower(data[,idlist])
  iddata_mysql <- iddata_mysql[,(!(colnames(iddata_mysql) %in% colnames(data))) | colnames(iddata_mysql) == "Compound ID"]
  
  predata <- merge(x = data,y = iddata_mysql,by.x = "newlist",by.y = "Compound ID",all.x = T,sort = F)[,-1,drop = F]
  
  return(predata)
}

#' @export
getmetamz <- function(data){
  library(dplyr)
  hmdbdata3 <- data
  
  for (i in 1:dim(hmdbdata3)[1]) {
    
    SMILES <- hmdbdata3[i,"smiles"]
    
    if(is.na(SMILES)){
      next
    }
    
    strsplit(SMILES,split = "+",fixed = T)
    strsplit(SMILES, "-", fixed =T)
    
    posnum <- strsplit(paste("a",SMILES,"a"), "+", fixed =T) %>% lapply(., length) %>% unlist()
    negnum <- strsplit(paste("a",SMILES,"a"), "-", fixed =T) %>% lapply(., length) %>% unlist()
    totalnum <- posnum-negnum
    if(totalnum == 0){
      hmdbdata3[i,"mode"] <- NA
      hmdbdata3[i,"charge"] <- NA
      hmdbdata3[i,"Adduct"] <- NA
    }else if(totalnum > 0){
      hmdbdata3[i,"mode"] <- "pos"
      hmdbdata3[i,"charge"] <- totalnum
      hmdbdata3[i,"Adduct"] <- "M+"
    }else{
      hmdbdata3[i,"mode"] <- "neg"
      hmdbdata3[i,"charge"] <- abs(totalnum)
      hmdbdata3[i,"Adduct"] <- "M-"
    }
  }
  
  hmdbdata4 <- hmdbdata3[is.na(hmdbdata3$Adduct),]
  hmdbdata5 <- hmdbdata3[!is.na(hmdbdata3$Adduct),]
  hmdbdata5 <- hmdbdata5[hmdbdata5$charge <= 1,]
  hmdbdata5[,"mz"] <-  hmdbdata5$Exact_Mass
  hmdbdata5[,"adductElement"] <- NA  
  hmdbdata5[,"nM"] <- NA
  
  # adduct <- getmysqldata(dbname = "meta",
  #                        table = "spatialadduct",
  #                        wherename = "usetype",
  #                        wheredata = 1)
  adduct <- getmysqldata(dbname = "meta",
                         table = "spatialadduct")
  adduct <- adduct[adduct$Type != "isotope",]
  
  for ( i in 1:dim(adduct)[1]) {
    hmdbdata6 <- hmdbdata4[!is.na(hmdbdata4$Mass),]
    if(dim(hmdbdata6)[1]==0){
      next
    }
    hmdbdata6$charge <- adduct[i,"charge"]
    hmdbdata6$mode <- adduct[i,"Type"]
    hmdbdata6$Adduct <- adduct[i,"Name"]
    hmdbdata6[,"mz"] <- (hmdbdata6$Mass*adduct[i,"nM"])/adduct[i,"charge"]+adduct[i,"Mass"]
    hmdbdata6[,"adductElement"] <- adduct[i,"Element"]
    hmdbdata6[,"nM"] <- adduct[i,"nM"]
    hmdbdata5 <- rbind(hmdbdata5,hmdbdata6)
  }
  
  for ( i in 1:dim(hmdbdata5)[1]) {
    if(is.na(hmdbdata5[i,"adductElement"])){
      hmdbdata5[i,"Formula-adduct"] <- hmdbdata5[i,"Formula"]
    }else{
      hmdbdata5[i,"Formula-adduct"] <- getnewformula(formula = hmdbdata5[i,"Formula"],
                                                     adductElement = hmdbdata5[i,"adductElement"],
                                                     nM = hmdbdata5[i,"nM"])
    }
  }
  hmdbdata5 <- hmdbdata5[!is.na(hmdbdata5[,"Formula-adduct"]),]
  
  return(hmdbdata5)
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",default = NULL, 
                      help = "表达数据矩阵文件,如果是xlsx文件以`数据矩阵.xlsx 数据矩阵`形式传参",
                      nargs = "+")
  parser$add_argument("-i","--idlist", default = "Compound ID", 
                      help="表格中id对应列名")
  parser$add_argument("-m","--metaid",default = NULL, 
                      help = "代谢物id编号",
                      nargs = "+")
  parser$add_argument("-df","--databasefrom",default =  c("inchikey","metlin","hmdb","lipidmaps","cas","kegg","英文名"), 
                      help = "数据查询数据库",
                      choices = c("inchikey","metlin","hmdb","lipidmaps","cas","kegg","英文名","中文名","来源"),
                      nargs = "+")
  parser$add_argument("-v","--vague", default = F, action = "store_true", 
                      help="是否模糊匹配")
  parser$add_argument("-d","--distinct", default = F, action = "store_true", 
                      help="是否去重")

  
  args <- parser$parse_args()
  if(is.null(args$filename)){
    args$data <- data.frame("Compound ID" = args$metaid,check.names = F)
    args$idlist <- "Compound ID"
  }else{
    args$data <- readdata(filename = filename)
  }
  args$filename <- NULL
  args$metaid <- NULL
  
  data <- do.call(what = getmetainfo, args = args)
  
  savexlsx(data = data,filename = "代谢物信息.xlsx",sheet = "代谢物信息")
  
}
