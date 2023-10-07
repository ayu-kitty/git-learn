#!/opt/conda/bin/Rscript

#' 获取空代最终定性结果
#'
#' @param savepeakpath peak信息路径
#' @param saveannopath 定性信息路径
#' @param moderange 正负离子模式
#'
#' @export
getfinalqualitative <- function(savepeakpath = "./sample/peak/",
                                saveannopath = "./sample/qualitative/",
                                moderange = c("neg","pos"),
                                deldup = T){
  
  qualitativedata2 <- NULL
  sheetname <- getsheetname(filename = paste0(savepeakpath,"/样本离子统计.xlsx"))
  
  for (mode in moderange) {
    if(!(mode %in% sheetname)){
      next
    }
    
    referencemz <- readdata(filename = paste0(savepeakpath,"/样本离子统计.xlsx"),sheet = mode)
    referencemz <- referencemz[referencemz$sampleselect,]
    referencemz <- referencemz$mz
    qualitativedata <- readdata(filename = paste0(saveannopath,"/adducts.xlsx"),sheet = mode)
    qualitativedata <- qualitativedata[qualitativedata$mz %in% referencemz,]
    qualitativedata <- qualitativedata[!(qualitativedata$polymer %in% "C2H4O"),]
    qualitativedata <- qualitativedata[!is.na(qualitativedata$score),]
    qualitativedata <- qualitativedata[is.na(qualitativedata$isotope),]
    qualitativedata <- qualitativedata[order(qualitativedata$intensity,decreasing = T),]
    qualitativedata <- qualitativedata[!duplicated(qualitativedata$formula),]
    qualitativedata <- qualitativedata[order(qualitativedata$mz),]
    qualitativedata[,"mode"] <- mode
    qualitativedata2 <- rbind(qualitativedata2,qualitativedata)
  }
  
  if(deldup){
    qualitativedata2 <- qualitativedata2[order(qualitativedata2$score,decreasing = T),]
    qualitativedata2 <- qualitativedata2[!duplicated(qualitativedata2$formula),]
    qualitativedata2 <- qualitativedata2[order(qualitativedata2$mz),]
  }
  
  returndata <- list()
  for (mode in moderange) {
    qualitativedata <- qualitativedata2[qualitativedata2$mode == mode,]
    if(dim(qualitativedata)[1] == 0){
      returndata[[mode]] <- NULL
      next
    }
    
    returndata[[mode]] <- qualitativedata$mz
    annodata <- readdata(paste0(saveannopath,"/anno.xlsx"),sheet =  mode)
    
    peakdata2 <- qualitativedata[,c("mz","adducts","formula","ppm")]
    colnames(peakdata2) <- c("mz","Adduct","Formula","ppm")
    peakdata2[,"ID"] <- 1:dim(peakdata2)[1]
    annodata3 <- NULL
    
    for (i in 1:dim(peakdata2)[1]) {
      if(is.na(peakdata2[i,"Formula"])){
        next
      }
      annodata2 <- annodata[annodata$Formula == peakdata2[i,"Formula"], ]
      if(dim(annodata2)[1]==0){
        next
      }
      
      annodata2[,"ID"] <- peakdata2[i,"ID"]
      annodata2[,"ID-num"] <- paste0(peakdata2[i,"ID"],"-",1:dim(annodata2)[1])
      annodata2[,"Adduct"] <- peakdata2[i,"Adduct"]
      annodata2[,"Formula"] <- peakdata2[i,"Formula"]
      annodata2[,"ppm"] <- peakdata2[i,"ppm"]
      annodata2[,"mz"] <- peakdata2[i,"mz"]
      annodata3 <- rbind(annodata3,annodata2)
      
      annodata2[is.na(annodata2)] <- ""
      # peakdata2[i, "Metabolites"] <- paste(annodata2$Metabolites, collapse = ";\r\n")
      # peakdata2[i, "Compound ID"] <- paste(annodata2$`Compound ID`, collapse = ";\r\n")
      # peakdata2[i, "Super Class"] <- paste(annodata2$`Super Class`, collapse = ";\r\n")
      # peakdata2[i, "Class"] <- paste(annodata2$`Class`, collapse = ";\r\n")
      # peakdata2[i, "Sub Class"] <- paste(annodata2$`Sub Class`, collapse = ";\r\n")
      # peakdata2[i, "KEGG"] <- paste(annodata2$`KEGG`, collapse = ";\r\n")
      peakdata2[i, "Metabolites"] <- paste(annodata2$Metabolites, collapse = "; ")
      peakdata2[i, "Compound ID"] <- paste(annodata2$`Compound ID`, collapse = "; ")
      peakdata2[i, "Super Class"] <- paste(annodata2$`Super Class`, collapse = "; ")
      peakdata2[i, "Class"] <- paste(annodata2$`Class`, collapse = "; ")
      peakdata2[i, "Sub Class"] <- paste(annodata2$`Sub Class`, collapse = "; ")
      peakdata2[i, "KEGG"] <- paste(annodata2$`KEGG`, collapse = "; ")
      peakdata2[i, "level"] <- annodata2$`level`[1]
      peakdata2[i, "num"] <- dim(annodata2)[1]
    }
    
    peakdata2 <- peakdata2[,c("ID","mz","Formula","Adduct","Metabolites","Compound ID","Super Class","Class","Sub Class","KEGG","level","num","ppm")]
    annodata3 <- annodata3[,c("ID","ID-num","Compound ID","Metabolites","Super Class","Class","Sub Class","KEGG","Formula","Adduct","mz","level","ppm")]
    peakdata2[,"mz"] <- format(peakdata2[,"mz"], nsmall = 5, trim = T)
    annodata3[,"mz"] <- format(annodata3[,"mz"], nsmall = 5, trim = T)
    
    savexlsx1(data = peakdata2,
              filename = paste0(saveannopath,"Qualitative.xlsx"),
              sheet = mode)
    savexlsx1(data = annodata3,
              filename = paste0(saveannopath,"Qualitative.xlsx"),
              sheet = paste0(mode,"-all"))
    
  }
  
  return(returndata)
}


#' 获取空代最终定性结果
#'
#' @param savepeakpath peak信息路径
#' @param saveannopath 定性信息路径
#' @param moderange 正负离子模式
#'
#' @export
getfinalqualitative_2023 <- function(saveannopath = "./sample/qualitative/",
                                     moderange = c("neg","pos"),
                                     deldup = T){
  
  qualitativedata2 <- NULL
  sheetname <- getsheetname(filename = paste0(saveannopath,"/adducts.xlsx"))
  
  for (mode in moderange) {
    if(!(mode %in% sheetname)){
      next
    }
    qualitativedata <- readdata(filename = paste0(saveannopath,"/adducts.xlsx"),sheet = mode)
    qualitativedata <- qualitativedata[qualitativedata$sampleselect,]
    qualitativedata[,"mode"] <- mode
    qualitativedata2 <- rbind(qualitativedata2,qualitativedata)
  }
  
  if(deldup){
    qualitativedata2 <- qualitativedata2[order(qualitativedata2$score,decreasing = T),]
    qualitativedata2 <- qualitativedata2[!duplicated(qualitativedata2$formula),]
    qualitativedata2 <- qualitativedata2[order(qualitativedata2$mz),]
  }
  
  returndata <- list()
  for (mode in moderange) {
    qualitativedata <- qualitativedata2[qualitativedata2$mode == mode,]
    if(dim(qualitativedata)[1] == 0){
      returndata[[mode]] <- NULL
      next
    }
    
    returndata[[mode]] <- qualitativedata$mz
    annodata <- readdata(paste0(saveannopath,"/anno.xlsx"),sheet =  mode)
    annodata <- getmetainfo(annodata,idlist = "cid",cid = T,needlist = c("CAS","HMDB","METLIN","Lipidmaps","KEGG","PubChem","ChEBI",
                                                                         "smiles","InChIKey","Mass","Formula",
                                                                         "Super Class","Class","Sub Class"))
    colnames(annodata)[colnames(annodata) == "mz"] <- "Adduct Mass"
    
    peakdata2 <- qualitativedata[,c("mz","adducts","formula","ppm")]
    colnames(peakdata2) <- c("mz","Adduct","Formula","ppm")
    peakdata2[,"ID"] <- 1:dim(peakdata2)[1]
    annodata3 <- NULL
    
    for (i in 1:dim(peakdata2)[1]) {
      if(is.na(peakdata2[i,"Formula"])){
        next
      }
      annodata2 <- annodata[annodata$Formula == peakdata2[i,"Formula"], ]
      if(dim(annodata2)[1]==0){
        next
      }
      
      annodata2[,"ID"] <- peakdata2[i,"ID"]
      annodata2[,"ID-num"] <- paste0(peakdata2[i,"ID"],"-",1:dim(annodata2)[1])
      # annodata2[,"Adduct"] <- peakdata2[i,"Adduct"]
      # annodata2[,"Formula"] <- peakdata2[i,"Formula"]
      annodata2[,"ppm"] <- peakdata2[i,"ppm"]
      annodata2[,"mz"] <- peakdata2[i,"mz"]
      annodata3 <- rbind(annodata3,annodata2)
      
      annodata2[is.na(annodata2)] <- ""
      # peakdata2[i, "Metabolites"] <- paste(annodata2$Metabolites, collapse = ";")
      # peakdata2[i, "Compound ID"] <- paste(annodata2$`Compound ID`, collapse = ";")
      # peakdata2[i, "Super Class"] <- paste(annodata2$`Super Class`, collapse = ";")
      # peakdata2[i, "Class"] <- paste(annodata2$`Class`, collapse = ";")
      # peakdata2[i, "Sub Class"] <- paste(annodata2$`Sub Class`, collapse = ";")
      # peakdata2[i, "KEGG"] <- paste(annodata2$`KEGG`, collapse = ";")
      for ( needlist in c("Metabolites","cid","CAS","HMDB","METLIN","Lipidmaps","KEGG","PubChem","ChEBI",
                          "smiles","InChIKey","Super Class","Class","Sub Class")) {
        peakdata2[i,needlist] <- paste(annodata2[,needlist], collapse = "; ")
      }
      
      peakdata2[i, "level"] <- annodata2$`level`[1]
      peakdata2[i, "Mass"] <- annodata2$`Mass`[1]
      peakdata2[i, "Adduct Mass"] <- annodata2$`Adduct Mass`[1]
      peakdata2[i, "num"] <- dim(annodata2)[1]
    }
    
    peakdata2 <- peakdata2[,c("ID","mz","Adduct Mass","ppm","Formula","Adduct","Metabolites","cid","HMDB","METLIN","Lipidmaps","KEGG","ChEBI","PubChem","CAS","smiles","InChIKey","Super Class","Class","Sub Class","level","num")]
    annodata3 <- annodata3[,c("ID","ID-num","mz","Adduct Mass","ppm","Metabolites","cid","HMDB","METLIN","Lipidmaps","KEGG","ChEBI","PubChem","CAS","smiles","InChIKey","Super Class","Class","Sub Class","Formula","Adduct","level")]
    peakdata2[,"mz"] <- format(peakdata2[,"mz"], nsmall = 5, trim = T)
    annodata3[,"mz"] <- format(annodata3[,"mz"], nsmall = 5, trim = T)
    
    savexlsx1(data = peakdata2,
              filename = paste0(saveannopath,"Qualitative.xlsx"),
              sheet = mode)
    savexlsx1(data = annodata3,
              filename = paste0(saveannopath,"Qualitative.xlsx"),
              sheet = paste0(mode,"-all"))
    
  }
  
  return(returndata)
}