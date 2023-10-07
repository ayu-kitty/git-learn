#' @export
lipid_predeal <- function(filename_neg=NA, filename_pos=NA,samplename=NULL){
  
  pacman::p_load(dplyr,tidyverse)
  
  print("开始预处理")
  neg <- predeal.lipid(filename=filename_neg,mode="neg")
  neg <- neg %>% dplyr::rename(setNames(samplename[,"negname"], samplename[,"samplename"]))
  pos <- predeal.lipid(filename=filename_pos,mode="pos")
  pos <- pos %>% dplyr::rename(setNames(samplename[,"posname"], samplename[,"samplename"]))
  data_merged <- rbind( pos,neg) %>%
    arrange(desc(`Grade-D`)) %>%
    arrange(desc(`Grade-C`)) %>%
    arrange(desc(`Grade-B`)) %>%
    arrange(desc(`Grade-A`))
  data_merged_unique <- data_merged[!duplicated(data_merged$LipidIon),]
  
  # 增加ID列，并对列名进行更改
  data_final <- mutate(data_merged_unique,
                       ID= 1:dim(data_merged_unique)[1]) %>% 
    plyr::rename(c("LipidIon"="Metabolites",
                   "Retention time"="Retention time (min)",
                   "IonFormula"="Formula",
                   "CalcMz"="m/z",
                   "Ionmode"="Ion mode",
                   "Class"="Sub Class",
                   "LipidGroup"="Adducts"))
  data_final[,"Total Carbon"] <- gsub(pattern = ":.*",replacement = "",x = data_final[,"Adducts"])
  data_final[,"Total Carbon"] <- as.numeric(gsub(pattern = ".*[a-zA-z(]",replacement = "",x = data_final[,"Total Carbon"]))
  data_final[,"Total Unsaturation"] <- gsub(pattern = ".*:",replacement = "",x = data_final[,"Adducts"])
  data_final[,"Total Unsaturation"] <- as.numeric(gsub(pattern = "[a-zA-z)].*",replacement = "",x = data_final[,"Total Unsaturation"]))
  data_final[,"Connection Method1"] <- gsub(pattern = ".*[(]",replacement = "",x = data_final[,"Adducts"])
  data_final[grepl(pattern = "e",data_final[,"Connection Method1"]),"Connection Method"] <- "e"
  data_final[grepl(pattern = "p",data_final[,"Connection Method1"]),"Connection Method"] <- "p"
  
  data_final[,"Adducts"] <- gsub(pattern = ".*[)]",replacement = "M",x = data_final[,"Adducts"])
  dataclass <- readdata(databasepath("database/database/qualitative/UtargetLipid/Lipid_match.xlsx"),sheet="ULipid")
  data_final <- merge(data_final,dataclass,by = "Sub Class",sort = F)
  
  data_final_info <- select(data_final,!samplename[,"samplename"])
  info <- getmetainfo(data = data_final[,"Metabolites",drop = F],idlist = "Metabolites",
                      needlist = c("cid","CAS","HMDB","METLIN","Lipidmaps","KEGG","PubChem","ChEBI","smiles","InChIKey"))
  data_final_info <- merge(data_final_info,info,by = "Metabolites",sort = F)
  data_final_info <- data_final_info[,c("ID","m/z","Retention time (min)","Ion mode","Metabolites",
                                        "cid","CAS","HMDB","METLIN","Lipidmaps","KEGG","PubChem","ChEBI","smiles","InChIKey",
                                        "Class","Sub Class","Total Carbon","Total Unsaturation","Connection Method",
                                        "FattyAcid","FA1","FA2","FA3","FA4","Adducts","Formula","Delta(ppm)","mScore","Grade")]
  data_final_data <-  select(data_final,c("ID",samplename[,"samplename"]))
  data_final <- merge(data_final_info,data_final_data,by = "ID")
  
  print("数据预处理完成")
  
  return(data_final)
}


#' @export
predeal.lipid <- function(filename,mode){
  suppressWarnings(library(stringr))
  
  if(is.na(filename)){
    data_final <-NULL
  }else{
    data2 <- read.csv(filename,header = F, sep = "\t", encoding = 'utf-8', check.names = F,stringsAsFactors = F)
    num <- which(data2[,1] == "Rej.")
    data <- read.csv(filename,header = T, sep = "\t", encoding = 'utf-8', check.names = F,stringsAsFactors = F,skip = num)
    num <- which(data2[,1] == "#normalize base:")
    dataname <- data2[1:num-1,,drop = F]	
    dataname <- tidyr::separate(data = dataname,sep = ":",col = 1,into = c("raw","name"))
    dataname$name <- gsub(pattern = ".raw$",replacement = "",x = dataname$name)
    dataname$raw <- gsub(pattern = "#",replacement = "Area",x = dataname$raw)
    
    data <- data[data[,"Rej."]!=1,]
    data <- dplyr::rename(data,"drop"=length(colnames(data)))
    
    for ( i in 1:dim(data)[1]) {
      a<-unlist(strsplit(data$IonFormula[i],split = " "))
      b1<-str_extract(string=a,pattern = "[a-zA-Z]+")
      b2<-as.numeric(str_extract(string=a,pattern = "[0-9]+"))
      
      #2022.11.2 "+H$"统一修正为"//+H$"
      #2022.11.15 "新增+CH3COO加和离子处理"
      if(length(grep(data$LipidIon[i],pattern = "\\+H$"))!=0){
        b2[b1=="H"]<- b2[b1=="H"]-1
      }else if(length(grep(data$LipidIon[i],pattern = "\\+Na$"))!=0){
        b2[b1=="Na"]<- b2[b1=="Na"]-1
      }else if(length(grep(data$LipidIon[i],pattern = "\\+NH4$"))!=0){
        b2[b1=="H"]<- b2[b1=="H"]-4
        b2[b1=="N"]<- b2[b1=="N"]-1
      }else if(length(grep(data$LipidIon[i],pattern = "-H$"))!=0){
        b2[b1=="H"]<- b2[b1=="H"]+1
      }else if(length(grep(data$LipidIon[i],pattern = "-2H$"))!=0){
        b2[b1=="H"]<- b2[b1=="H"]+2
      }else if(length(grep(data$LipidIon[i],pattern = "\\+HCOO$"))!=0){
        b2[b1=="H"]<- b2[b1=="H"]-1
        b2[b1=="C"]<- b2[b1=="C"]-1
        b2[b1=="O"]<- b2[b1=="O"]-2
      }else if(length(grep(data$LipidIon[i],pattern = "\\+CH3COO$"))!=0){
        b2[b1=="H"]<- b2[b1=="H"]-3
        b2[b1=="C"]<- b2[b1=="C"]-2
        b2[b1=="O"]<- b2[b1=="O"]-2
      }
      
      data[i,"IonFormula"]<-paste(b1,b2,sep ="",collapse =" " )
    }
    
    deltappm<-grep("^Delta\\(ppm\\)*",colnames(data),value = T)
    data[,"Delta(ppm)"]<-apply(data[,deltappm], 1, function(x){
      setdiff(sort(unique(as.numeric(as.matrix(x))),decreasing = TRUE),1000000)[1]
    })
    mScore<-grep("mScore*",colnames(data),value = T)
    data[,"mScore"]<-apply(data[,mScore], 1, max)
    grade<-grep("Grade*",colnames(data),value = T)
    data[,"Grade"]<-apply(data[,grade],1,function(x)sort(setdiff(unique(x),""))[1])
    vari.name1 <- c('LipidIon', 'LipidGroup', 'Class', 'FattyAcid', 'FA1', 'FA2', 'FA3', 'FA4',
                    'CalcMz', 'IonFormula',"Delta(ppm)","mScore","Grade")
    if ('FA4' %in% colnames(data)){
      data1 <- select(data, all_of(vari.name1))
    }else if('FA3' %in% colnames(data)){
      data1 <- mutate(data, FA4='') %>%　select(all_of(vari.name1))
    }else if('FA2' %in% colnames(data)){
      data1 <- mutate(data, FA3='',FA4='') %>%　select(all_of(vari.name1))
    }else{
      data1 <- mutate(data, FA2='',FA3='',FA4='') %>%　select(all_of(vari.name1))
    }
    
    if (mode=="neg"){
      data1$LipidIon <- gsub('-H', '', data1$LipidIon) %>%
        gsub(pattern = '-2H', replacement = '') %>%
        gsub(pattern = "+HCOO", replacement = '', fixed = T)%>%
        gsub(pattern = "+CH3COO", replacement = '', fixed = T)
      data1[,"Ionmode"] <- "neg"
    }else{
      data1$LipidIon <- gsub('+H', '', data1$LipidIon, fixed = T) %>%
        gsub(pattern = '+Na', replacement = '', fixed = T) %>%
        gsub(pattern = "+NH4", replacement = '', fixed = T)
      data1[,"Ionmode"] <- "pos"
    }
    # vari.name2 <- c('Area[', 'Rt[', 'Grade[')
    
    # 归一化
    func <- function(x, na.rm = FALSE)(10000*x/sum(x, na.rm = FALSE))
    data2 <- select(data, starts_with('Area[')) %>% transmute_all(.funs = func)
    for ( i in 1:dim(dataname)[1]) {
      colnames(data2)[colnames(data2)==dataname[i,"raw"]] <-dataname[i,"name"]
    }
    
    # RT求平均值
    data3 <- data.frame("Retention time" = select(data, starts_with('Rt[')) %>% rowMeans(), check.names = F)
    
    # 按Grade排序
    mat <- function(z){
      sumz<-c(
        "Grade-A"=length(z[z=='A']),
        "Grade-B"=length(z[z=='B']),
        "Grade-C"=length(z[z=='C']),
        "Grade-D"=length(z[z=='D'])
      )
    }
    data4 <- t(apply(select(data, starts_with('Grade[')),MARGIN = 1, FUN = mat))
    
    
    data_bind <- cbind(data1, data3, data4,data2) %>%
      arrange(desc(`Grade-D`)) %>%
      arrange(desc(`Grade-C`)) %>%
      arrange(desc(`Grade-B`)) %>%
      arrange(desc(`Grade-A`))
    
    # 去重LipidIon
    data_final <-  data_bind[!duplicated(data_bind$LipidIon),]
    #RMETA2::savexlsx1(data =data_final,name = "数据矩阵.xlsx",sheet = mode)
  }
  return(data_final)
}



