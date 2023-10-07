#! /opt/conda/bin/Rscript

#' 1.get raw data 
#' @param tlcms rawData
#' @export
getRawData_target <- function(tlcms=NA,rawData=NA){
  library(openxlsx)
  options(scipen = 200)
  data <- list()
  #可优化（简化）
  if(is.na(tlcms)){
    tlcms <- list.files(pattern = "TLCMS")
  }else{
    tlcms=tlcms
  }
  sheets <- getSheetNames(tlcms)
  for (sheet in sheets) {
    if(grepl("基本参数",sheet)){
      info <- openxlsx::read.xlsx(tlcms,sheet)
      data[["info_base"]] <- info
      proName <- paste0(
        info[info[1]=="项目编号",2],"-",
        info[info[1]=="客户单位",2],"-",
        info[info[1]=="物种",2],"-",
        info[info[1]=="样本类型",2],"-",
        info[info[1]=="靶向类型",2],
        "（项目报告）"
      )
      data[["proName"]]=proName
    }else if(grepl("样本信息",sheet)){
      info <- openxlsx::read.xlsx(tlcms,sheet)
      data[["info"]] <- info
    }else if(grepl("绝对定量",sheet)){
      info <- openxlsx::read.xlsx(tlcms,sheet)
      v1=1.00
      n=1.00
      v2=1.00
      
      if(dim(info)[2]==4){
        v1=as.double(info[3])
        n=as.double(info[4])
      }else if(dim(info)[2]==5){
        v1=as.double(info[3])
        n=as.double(info[4])
        v2=as.double(info[5])
      }
      info_cal <- list(v1,n,v2)
      names(info_cal) <- c("V1","N","V2")
      data[["info_cal"]] <- info_cal
      
    }else if(grepl("比较组信息",sheet)){
      info <- openxlsx::read.xlsx(tlcms,sheet)
      info1 <- as.vector(unlist(info[3]))
      if(!is.na(info1[1])){
        data[["compare_list"]] <- info1
		data[["sample_group"]] <- info
      }else{
        data[["compare_list"]] <- NA
		data[["sample_group"]] <- NA
      }
    }else{
      data=data
    }
  }
  
  # get raw data
  data_raw=NA
  if(endsWith(rawData,"txt")){
    data_raw <- read.csv(file = rawData,
                         check.names = "F",
                         sep = "\t",
                         encoding = "UTF-8",
                         stringsAsFactors = F,
                         na.strings = c("< 0","N/A","<2 points","degenerate"))
  }else if(endsWith(rawData,"csv")){
    data_raw <- read.csv(file = rawData,
                         check.names = "F",
                         sep = ",",
                         encoding = "UTF-8",
                         stringsAsFactors = F,
                         na.strings = c("< 0","N/A","N/F","<2 points","degenerate"))
  }else if(endsWith(rawData,"xlsx")){
    
    data_raw <- read.xlsx(xlsxFile=rawData,
                          sheet = 1,
                          sep.names = " ",
                          na.strings = c("< 0","N/A","<2 points","degenerate"))
    if(length(getSheetNames(rawData))>1){
      for (sheet in getSheetNames(rawData)[-1]) {
        temp1 <- read.xlsx(xlsxFile=rawData,
                           sheet=sheet,
                           sep.names = " ",
                           na.strings = c("< 0","N/A","<2 points","degenerate"))
        data_raw <- rbind(data_raw,temp1)
      }
    }else{
      data_raw=data_raw
    }

  }else{
      print("当前文件格式不支持")
  }
  data_raw[is.na(data_raw)] <- 0
  data[["raw"]] <- data_raw
  return(data)
}

##'2 processing raw data
#' @export
rawDataProcess_target <- function (data = data,
                            info_main=NA,
                            info_del=NA,
                            info_filter=NA) {
  data_proccessed = NA
  cols_del=seq(1,length(info_del),2)
  cols_filter=seq(1,length(info_filter),2)
  #选取非重复列
  all_cols <- unique(c(info_main,info_del[cols_del],info_filter[cols_filter]))
  all_cols <- all_cols[!is.na(all_cols)]
  data_proccessed <- data$raw[all_cols]
  #删除
  if(is.na(info_del[1])){
    data_proccessed <- data_proccessed
  }else{
    #取反
    for (del1 in cols_del) {
      data_proccessed <- data_proccessed[which(data_proccessed[info_del[del1]]!=info_del[del1+1]),]
    }
  }
  #筛选
  if(is.na(info_filter[1])){
    data_proccessed = data_proccessed
  }else{
    temp1 <- NA
    for (filter1 in cols_filter) {
      if(filter1==1){
        # filter1=1 #测试用语句
        temp1 <- data_proccessed[which(data_proccessed[info_filter[filter1]]==info_filter[filter1+1]),]
      }else{
        # filter1=3 #测试用语句
        temp1 <- data_proccessed[which(data_proccessed[info_filter[filter1]]!=as.numeric(info_filter[filter1+1])),]
        }
      }
    temp1 <- unique(temp1[info_main[2]])
    data_proccessed <- merge(temp1,data_proccessed,
                             by = info_main[2],
                             all.x = T)
    }
  data_proccessed <- data_proccessed[info_main]
  colnames(data_proccessed) <- c("Sample Name", "Metabolites","Area", "Concentration")
  for (i in 1:(dim(data_proccessed)[1])) {
    meta_old <- data_proccessed[i,2]
    if(grepl("α",meta_old)){
      meta_new <- gsub("α","Alpha",meta_old)
      data_proccessed[i,2] <- meta_new
    }else if(grepl("β",meta_old)){
      meta_new <- gsub("β","Beta",meta_old)
      data_proccessed[i,2] <- meta_new
    }else if(grepl("γ",meta_old)){
      meta_new <- gsub("γ","Gamma",meta_old)
      data_proccessed[i,2] <- meta_new
    }else{
      data_proccessed[i,2] <- meta_old
    }
  }
  data[["raw_processed"]] <- data_proccessed
  return(data)
}
##'3.get area and cal data from procceded file
#' @export
getAreaCal_target <- function(data = data,qc_name="STD-QC",drop_tails="_[1-9]$") {
  # data=data
  # qc_name="STD-QC"
  library(tidyr)
  library(dplyr)
  
  library(data.table)
  options(scipen = 200)
  
  names_old <- c("Metabolites",data$info[which(!is.na(data$info[4])),1])
  names_new <- c("Metabolites",data$info[which(!is.na(data$info[4])),2])
  qc_old <- NA
  qc_new <- NA
  if(dim(data$info)[2]==4){
    qc_old <- c("Metabolites",data$info[which(data$info[3]==qc_name),1])
    qc_new <- c("Metabolites",data$info[which(data$info[4]==qc_name),2])
  }else{
    qc_old <- NA
    qc_new <- NA
  }
  
  area <- data$raw_processed%>%
    pivot_wider(names_from = "Sample Name",
                 id_cols = "Metabolites",
                values_from = "Area")%>%
    select(all_of(names_old))
  
  cal <- data$raw_processed%>%
    pivot_wider(names_from = "Sample Name",
                id_cols = "Metabolites",
                values_from = "Concentration")%>%
    select(all_of(names_old))
  
  meta_sub <- function(data1){
    return(gsub(drop_tails,"",data1))
  }
  # drop_tails=NA
  if(!is.na(drop_tails)){
    area["Metabolites"] <- apply(area["Metabolites"],1,meta_sub)
    cal["Metabolites"] <- apply(cal["Metabolites"],1,meta_sub)
  }
  
  if(!any(is.na(qc_old))){
    qc <- data$raw_processed%>%
      pivot_wider(names_from = "Sample Name",
                  id_cols = "Metabolites",
                  values_from = "Area")%>%
      select(all_of(qc_old))
    qc["Metabolites"] <- apply(qc["Metabolites"],1,meta_sub)
  }else{
    qc <- NA
  }
  
  #代谢物名称修改
  setnames(area,names_old,names_new)
  setnames(cal,names_old,names_new)
  
  data[["qc"]] <- qc
  area[is.na(area)] <- 0
  data[["area"]] <- area
  cal[is.na(cal)] <- 0
  data[["cal"]] <- cal
  return(data)
}

##'4.get stable info of data
#' @export
getStable_target <- function(data = data,drop_tails=T){
  options(scipen = 200)
  library(dplyr)
  library(tidyr)
  library(data.table)
  temp <- data[["qc"]]
  
  temp1 <- rowMeans(temp[-1])
  temp2 <- apply(temp[-1], 1,sd)
  temp3 <- temp2/temp1*100
  stable <- data.frame(temp1,temp3)
  names(stable) <- c("Ave(STD)","Rsd(STD)")
  stable <- cbind(temp[1],stable)
  data[["stable"]] <- stable
  return(data)
}

##'5.calculate Con
#' @export
calCon_target <- function(data=data){
  options(scipen = 200)
  library(openxlsx)
  groups <- unlist(unique(data$info[which(!is.na(data$info[4])),3]))
  temp <- data$cal
  con_groups <- list()
  #group_L = c("ID","Metabolites")
  for (group in groups) {
    cal_cols <- data$info[data$info[3] == group,2]
	#group_L <- c(group_L,data$info[data$info[3] == group,3])
    cols_weight <- data$info[data$info[3] == group,dim(data$info)[2]]

    #计算tapply
    temp1 <- temp[cal_cols]*data$info_cal$V1*data$info_cal$N
    temp1 <- t(t(temp1)/cols_weight)
    con_groups <- c(con_groups,
                    list(temp1))
  }
  con <- data.frame(con_groups,check.names = F)
  #添加ID列,索引列
  #con <- cbind(ID=row.names(con),temp[1],con)
  con <- cbind(temp[1],con)
  # 将group_L转换为数据框，并使用con的列名
  #group_L_df <- data.frame(t(group_L))
  #colnames(group_L_df) <- colnames(con)
  #con <- rbind(group_L_df,con)
  con[is.na(con)] <- 0
  data[["con"]] <- con
  
  groups[is.na(groups)] <- 0
  names(con_groups) <- groups
  data[["con_groups"]] <- con_groups
  return(data)
}

##'save data by names of data
#' @export
saveData_target <- function(data=data,saveRsd=T){
  data1 <- list()
  data_names <- list()

  if(any(names(data)=="info") && !is.null(data[["info"]])){
    data1 <- c(data1,list(data$info[,2:3]))
    data_names <- c(data_names,"样本信息")
  }
    #area
  if(any(names(data)=="area") && !is.null(data[["area"]])){
    data1 <- c(data1,list(data[["area"]]))
    data_names <- c(data_names,"峰面积信息")
  }
    #stable
  if(saveRsd){
    data1 <- c(data1,list(data$stable))
    data_names <- c(data_names,"稳定性信息")
  }else{
    #do nothing
    data1 <- data1
  }
    #cal
  # if(any(names(data)=="cal") && !is.na(data[["cal"]])){
  if(any(names(data)=="cal") && !is.null(data[["cal"]])){
    data1 <- c(data1,list(data$cal))
    data_names <- c(data_names,"浓度信息")
  }
    #con
  # if(any(names(data)=="con") && !is.na(data[["con"]])){
  if(any(names(data)=="con") && !is.null(data[["con"]])){
    data1 <- c(data1,list(data[["con"]]))
    data_names <- c(data_names,"定量信息")
  }
    #比较组信息
  if(any(names(data)=="sample_group") && !is.null(data[["sample_group"]])){
    data1 <- c(data1,list(data[["sample_group"]]))
    data_names <- c(data_names,"比较组信息")
  }
    #save
  names(data1) <- data_names
  #openxlsx::write.xlsx(data1,"data.xlsx")
  #save(data1,file = "raw.RData")
  return(data1)
}

#' 从数据库根据Metabolites 获得对应的注释信息
#' @param predata 要注释的数据
#' @export
get_annotation <- function(
              predata = predata,
              protype = protype
  ){
  Metabolites_pubchem <- openxlsx::read.xlsx(paste0("/data/hstore2/database/database/靶向分析/",protype,"_cid.xlsx"), sheet = 1)
  #Metabolites_pubchem <- openxlsx::read.xlsx(paste0("/data/hstore4/database/靶向分析/",protype,"_cid.xlsx"), sheet = 1)
  #Metabolites和Name不区分大小写，便于匹配
  #predata$merge_key <- tolower(predata$Metabolites)
  #Metabolites_pubchem$merge_key <- tolower(iddata_mysql_info$Name)
  iddata_mysql <- merge(predata,Metabolites_pubchem,by = "Metabolites", all.x = T) 
  iddata_mysql_info <- getmysqldata(dbname = "cosa",
                                    table = "compound_identifier",
                                    tablelist = "cid,formula as Formula,pubchemid as PubChem,casNumber as CAS,hmdbid as HMDB,metlinid as METLIN,lipidmapsid as Lipidmaps,keggid as KEGG")
  iddata_mysql <- merge(iddata_mysql,iddata_mysql_info,by = "cid", all.x = T)

  iddata_mysql_inchikey <- getmysqldata(dbname = "cosa",
                                        table = "compound_structure",
                                        tablelist = "cid,smiles,inchikey as InChIKey",
                                        wherename = "cid",
                                        wheredata = iddata_mysql$cid)
  iddata_mysql <- merge(iddata_mysql,iddata_mysql_inchikey,by = "cid",all.x = T)
  
  iddata_mysql_chebiid <- getmysqldata(dbname = "cosa",
                                        table = "extra_hmdb",
                                        tablelist = "cid,chebiid as ChEBI",
                                        wherename = "cid",
                                        wheredata = iddata_mysql$cid)
  iddata_mysql <- merge(iddata_mysql,iddata_mysql_chebiid,by = "cid",all.x = T)  
  
  iddata_mysql_classification <- getmysqldata(dbname = "cosa",
                                              table = "compound_classification",
                                              tablelist = "cid,superclass as 'Super Class',class as Class,subclass as 'Sub Class'",
                                              wherename = "cid",
                                              wheredata = iddata_mysql$cid)

  iddata_mysql <- merge(iddata_mysql,iddata_mysql_classification,by = "cid",all.x = T)
  return (iddata_mysql)
}

##' 主流程
#' @export
datafech_targettlcms <-function(inputpath = "./",
                                tlcms = "./TLCMS-胆汁酸.xlsx",
								qcinfo = "./胆汁酸-标品信息表.xlsx",
								rawData = NULL

){
  library(dplyr)
  file_name <- basename(qcinfo)
  split_names <- strsplit(file_name, "-")[[1]]
  protype <- split_names[1]
  setwd(inputpath)
  #标品名称对照
  name <- openxlsx::read.xlsx(qcinfo, sheet = 1)
  data <- getRawData_target(tlcms = tlcms,rawData = rawData)

  #raw 中Component Name转换成英文名
  if (protype %in% c("胆汁酸")){
    name$鹿明内部编号 <- gsub("-", "_", name$鹿明内部编号)
    for(i in 1:nrow(name)){
      data$raw$`Component Name`[grepl(name$鹿明内部编号[i], data$raw$`Component Name`)] <- name$英文名[i]
    }
  }else if (protype %in% c("氨基酸")){
    data$raw$`Component Name` <- gsub("_1", "", data$raw$`Component Name`)
    name$英文名 <- gsub("L-","",name$英文名)
    name$英文名 <- gsub("4-Hydroxyproline","L-4-Hydroxyproline",name$英文名)
  }else if (protype %in% c("糖醇","TMAO","动物激素","短链脂肪酸","黄酮酚类")){
    data$raw$`Component Name` <- gsub("_1", "", data$raw$`Component Name`)
	data$raw$`Component Name` <- gsub("_3", "", data$raw$`Component Name`)
  }else if (protype %in% c("胆固醇","核苷")){
    for(i in 1:nrow(name)){
      data$raw$`Component Name`[grepl(name$鹿明内部编号[i], data$raw$`Component Name`)] <- name$英文名[i]
    }
  }else if (protype %in% c("有机酸")){
    data$raw$`Component Name` <- gsub("_1", "", data$raw$`Component Name`)
	data$raw$`Component Name` <- gsub("_2", "", data$raw$`Component Name`)
	data$raw$`Component Name` <- gsub("DL---phenyllactic acid", "DL-β--phenyllactic acid", data$raw$`Component Name`)
  }else if (protype %in% c("中长链脂肪酸")){
    data$raw$`Component Name` <- gsub("_1", "", data$raw$`Component Name`)
	data$raw$`Component Name` <- gsub("_2", "", data$raw$`Component Name`)
  }else if (protype %in% c("植物激素")){
    data$raw$`Component Name` <- gsub("-1", "", data$raw$`Component Name`)
	data$raw$`Component Name` <- gsub("-2", "", data$raw$`Component Name`)
    for(i in 1:nrow(name)){
      data$raw$`Component Name` <- ifelse(data$raw$`Component Name` == name$缩写[i], name$英文名[i], data$raw$`Component Name`)
    }
  }else if (protype %in% c("神经递质")){
     data <- data
  }else {
    stop(paste0(protype,"不在靶向分析类型中"))
  }

  #删除所有带IS后缀的行（例：LM-SCFA_IS02_2）
  data$raw <- data$raw[!grepl("_IS",data$raw$`Component Name`),]
  #实验部门不能统一规则，因此放宽限制，手动填写规则
  #info_main,前四列为样本名称、代谢物名称、峰面积值、计算值
  #info_del,用于删除的列与相应值
  #info_filter,用于筛选的列与相应值,其中第二个条件取反
  data <- rawDataProcess_target(data = data,
                       info_main = c("Sample Name","Component Name",
                                     "Area","Calculated Concentration"),
                       # info_del = NA,
                       info_del = c("Sample Name","blk ",
                                    "Sample Name","blk"),
                       info_filter = c("Sample Type","Standard",
                                       "Actual Concentration",0)
                       # info_filter = NA
  )

  #qc_name根据qc的名称取值
  #drop_tials去除相应的尾部后缀，通常后缀应为“_[1~3]$”格式(以“_”+“数字”为结尾)，若以"数字-"开头，drop_tials="^[1~3]-"
  data <- getAreaCal_target(data = data,drop_tails ="_[1~3]$")
  #获取稳定性信息
  data <- getStable_target(data = data)
  #计算含量结果
  data <- calCon_target(data = data)
  data <- saveData_target(data = data)

  #添加注释
  #表头顺序,去除QC
  samples <- data$样本信息$样本分析名[!grepl("QC",data$样本信息$样本分析名)]
  if (protype %in% c("胆汁酸","胆固醇")){
    #name中的β替换为beta
    name$英文名 <- gsub("α", "Alpha", name$英文名)
    name$英文名 <- gsub("β", "Beta", name$英文名)
    name$英文名 <- gsub("γ", "Gamma", name$英文名)
    anno_data1 <- merge(data$定量信息, name, by.x= "Metabolites", by.y="英文名", all.x=T) %>%
                 select(-鹿明内部编号,-化学式,-"分类（英文）",-分类,-Q1,-Q3,-RT) %>%
                 rename("Ion mode" = "离子模式")

    anno_data2 <- get_annotation(predata = anno_data1,protype = protype) %>%
                  select("Ion mode",Metabolites,cid,HMDB,METLIN,Lipidmaps,KEGG,ChEBI,PubChem,CAS,smiles,InChIKey,
				"Super Class",Class,"Sub Class",Formula,中文名,缩写,分子量,厂家,货号,纯度,all_of(samples))

  }else if (protype %in% c("黄酮酚类","中长链脂肪酸")){
    anno_data1 <- merge(data$定量信息, name, by.x= "Metabolites", by.y="英文名", all.x=T) %>%
                 select(-化学式,-"分类（英文）",-分类)
    anno_data2 <- get_annotation(predata = anno_data1,protype = protype) %>%
                  select(Metabolites,cid,HMDB,METLIN,Lipidmaps,KEGG,ChEBI,PubChem,CAS,smiles,InChIKey,
				"Super Class",Class,"Sub Class",Formula,中文名,分子量,厂家,货号,纯度,all_of(samples))
  }else if (protype %in% c("氨基酸","TMAO")){
    anno_data1 <- merge(data$定量信息, name, by.x= "Metabolites", by.y="英文名", all.x=T) %>%
                  select(-分子式)
    anno_data2 <- get_annotation(predata = anno_data1,protype = protype) %>%
                  select(Metabolites,cid,HMDB,METLIN,Lipidmaps,KEGG,ChEBI,PubChem,CAS,smiles,InChIKey,
				"Super Class",Class,"Sub Class",Formula,中文名,缩写,分子量,厂家,货号,纯度,all_of(samples))
  }else if(protype %in% c("糖醇")){
    anno_data1 <- merge(data$定量信息, name, by.x= "Metabolites", by.y="英文名", all.x=T) %>%
                  select(-分子式,RT) %>% rename("Ion mode" = "离子模式")
    anno_data2 <- get_annotation(predata = anno_data1,protype = protype) %>%
                  select("Ion mode",Metabolites,cid,HMDB,METLIN,Lipidmaps,KEGG,ChEBI,PubChem,CAS,smiles,InChIKey,
				"Super Class",Class,"Sub Class",Formula,中文名,分子量,厂家,货号,纯度,"特征碎片离子（m/z）",all_of(samples))
  }else if(protype %in% c("动物激素","核苷")){
    anno_data1 <- merge(data$定量信息, name, by.x= "Metabolites", by.y="英文名", all.x=T) %>%
                  select(-化学式,-分类)
    anno_data2 <- get_annotation(predata = anno_data1,protype = protype) %>%
                  select(Metabolites,cid,HMDB,METLIN,Lipidmaps,KEGG,ChEBI,PubChem,CAS,smiles,InChIKey,
				"Super Class",Class,"Sub Class",Formula,中文名,分子量,all_of(samples))
  }else if(protype %in% c("短链脂肪酸","植物激素")){
    anno_data1 <- merge(data$定量信息, name, by.x= "Metabolites", by.y="英文名", all.x=T)
    anno_data2 <- get_annotation(predata = anno_data1,protype = protype) %>%
                  select(Metabolites,cid,HMDB,METLIN,Lipidmaps,KEGG,ChEBI,PubChem,CAS,smiles,InChIKey,
				"Super Class",Class,"Sub Class",Formula,中文名,分子量,厂家,货号,纯度,all_of(samples))
  }else if(protype %in% c("神经递质")){
    anno_data1 <- merge(data$定量信息, name, by.x= "Metabolites", by.y="英文名", all.x=T) %>%
                  select(-CAS)
    anno_data2 <- get_annotation(predata = anno_data1,protype = protype) %>%
                  select(Metabolites,cid,HMDB,METLIN,Lipidmaps,KEGG,ChEBI,PubChem,CAS,smiles,InChIKey,
				"Super Class",Class,"Sub Class",Formula,中文名,缩写,all_of(samples))
  }else if(protype %in% c("有机酸")){
    anno_data1 <- merge(data$定量信息, name, by.x= "Metabolites", by.y="英文名", all.x=T) %>%
                  select(-分子式)
    anno_data2 <- get_annotation(predata = anno_data1,protype = protype) %>%
                  select(Metabolites,cid,HMDB,METLIN,Lipidmaps,KEGG,ChEBI,PubChem,CAS,smiles,InChIKey,
				"Super Class",Class,"Sub Class",Formula,中文名,分子量,all_of(samples))
  }else{
    stop(paste0(protype,"不在靶向分析类型中"))
  }
  
  anno_data2 <- cbind(ID=row.names(anno_data2),anno_data2) %>% 
               mutate(`Super Class` = if_else(is.na(`Super Class`), 'Unclassified', `Super Class`),
                      Class = if_else(is.na(Class), 'Unclassified', Class),
					  `Sub Class` = if_else(is.na(`Sub Class`), 'Unclassified', `Sub Class`))
  data <- c(data,list("注释矩阵" = anno_data2))
  openxlsx::write.xlsx(data,"data.xlsx")
  return (anno_data2)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  parser$add_argument("-ip","--inputpath", default = "./", help = "分析数据输入路径")
  parser$add_argument("-tl","--tlcms", default = "./TLCMS-胆汁酸.xlsx", help = "基本信息表", required = T)
  parser$add_argument("-qi","--qcinfo", default = "./胆汁酸-标品信息表.xlsx", help = "标准信息表")
  parser$add_argument("-d","--rawData", default = NULL, help = "下机数据", required = T)
  args <- parser$parse_args()
  result <- do.call(what = datafech_targettlcms,args = args)
}