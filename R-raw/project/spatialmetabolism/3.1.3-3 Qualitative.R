#!/opt/conda/bin/Rscript

#' 空代初步定性
#'
#' @param sample 选择样本名
#' @param areapath 选区路径
#' @param moderange 正负离子模式
#' @param maxnum 最大样本数
#' @param seed 随机选择种子
#' @param imzmlpath imzml数据路径
#' @param delsample 不运行的样本名
#' @param ... 见[GetimzmlpreQualitative()]
#'
#' @export
imzmlpreQualitative <- function(sample = NULL,
                                imzmlpath = "./sample/imzml-pre/",
                                areapath = "./sample/select/",
                                moderange = c("neg", "pos"),
                                maxnum = 10,
                                seed = 1111,
                                delsample = c("bg_data","qc_data"),
                                ...) {
  # 正负离子循环
  for (mode in moderange) {
    
    filename <- list.files(path = paste0(areapath, mode),
                           pattern = ".rds$",
                           full.names = F,
                           recursive = T)
    # 获取样品名
    samplename2 <- gsub(pattern = ".rds",
                        replacement = "",
                        x = filename)
    samplename <- gsub(pattern = ".*\\/", replacement = "", x = samplename2)
    
    slidename <- gsub(pattern = "\\/.*", replacement = "", x = samplename2)
    realslidename <- list.files(path = imzmlpath,pattern = paste0("-",mode,".imzML$"))
    realslidename <- gsub(pattern = paste0("-",mode,".imzML"),
                          replacement = "",
                          x = realslidename)
    samplename2 <- samplename2[slidename %in% realslidename]
    samplename <- samplename[slidename %in% realslidename]
    
    # 删除背景选区及QC
    if(!is.null(delsample)){
      for(delsample2 in delsample) {
        samplename2 <- samplename2[!grepl(pattern = paste0("^",delsample2,".*"),x = samplename)]
        samplename <- samplename[!grepl(pattern = paste0("^",delsample2,".*"),x = samplename)]
      }
    }
    
    if (any(duplicated(samplename))) {
      stop("样本名有重复，请核查")
    }
    
    if(!is.null(sample)){
      samplename2 <- samplename2[samplename %in% sample]
      samplename <- samplename[samplename %in% sample]
    }else if(length(samplename2) > maxnum){
      set.seed(seed =  seed)
      selectnum <- sample(1:length(samplename2),size = maxnum)
      samplename2 <- samplename2[selectnum]
      samplename <- samplename[selectnum]
    }
    
    if(length(samplename2) > 0){
      
      slidename <- gsub(pattern = "\\/.*", replacement = "", x = samplename2)
      filename <- paste0(imzmlpath,"/",slidename,"-",mode,".imzML")
      areards <- paste0(areapath,"/",mode,"/",samplename2,".rds")
      
      GetimzmlpreQualitative(filename = filename,
                             mode = mode,
                             areards = areards,
                             ...)
      
      gc(reset = TRUE)
      
    }else{
      print(paste0("在",areapath,"目录下未找到", mode, "模式的样本文件"))
    }
    
  }
}

#' 单一的imzml定性
#' 
#' @export
GetimzmlpreQualitative <- function(filename = NULL,
                                   areards = NULL,
                                   mode = "neg",
                                   # 读取参数
                                   mass.range = NULL,
                                   resolution = 5,
                                   units = "ppm",
                                   # 外置数据库
                                   outdatabse = NULL,
                                   # 对齐参数
                                   tolerance = 5,
                                   savepeakpath = "./sample/peak/",
                                   # 定性参数
                                   ppm = 5,
                                   polymernum = 3,
                                   polymerminmz = 300,
                                   # 保存参数
                                   saveannopath = "./sample/qualitative/",
                                   species_tissue = "小鼠",
                                   asp = 1,
                                   waters = F,
                                   ...){
  
  # setwd("/data/hstore2/test/2022-12-19空间代谢组数据测试/DLM202210597/")
  # roxygen2::roxygenize("/home/lujw/lmbio/lmbio-r")
  # 
  # filename = "sample/imzml/Lung-pos.imzML"
  # areards = "sample/select/pos/Lung/Lung.rds"
  # mode = "pos"
  # mass.range = c(70,1000)
  # resolution = 5
  # units = "ppm"
  # outdatabse = NULL
  # tolerance = 5
  # savepeakpath = "./sample/peak/"
  # ppm = 5
  # polymernum = 4
  # polymerminmz = 400
  # saveannopath = "./sample/qualitative/"
  # lcqualitativefrom = "./sample/qualitative/多项目非靶定性统计.xlsx"
  # manqualitativefrom = if(length(list.files(path = "./sample/qualitative/",pattern = "^人工数据库",full.names = T)) == 0){"人工数据库.xlsx"
  # }else{list.files(path = "./sample/qualitative/",pattern = "^人工数据库",full.names = T)}
  
  suppressMessages(library("Cardinal"))
  suppressMessages(library("mass2adduct"))
  suppressMessages(library("dplyr"))
  
  # 数据库获取
  if(!is.null(outdatabse)){
    # [x] 2023-03-29 添加中药数据库定性
    if(outdatabse == "tcm" | outdatabse == "tcme"){
      hmdbdata1 <- getmysqldata(dbname = "cosa",
                                table = "herb_annotation",
                                wherename = "is_animal_source",
                                wheredata = 0)
      hmdbdata1 <- hmdbdata1[,c("inchikey","compound_name","formula","exact_mass","smiles")]
      colnames(hmdbdata1) <- c("Compound ID","Metabolites","Formula","Exact_Mass","smiles")
      
      hmdbdata2 <- getmysqldata(dbname = "cosa",
                                table = "herb_ac")
      hmdbdata2 <- hmdbdata2[,c("Ingredient_inchikey","best_name","Ingredient_formula","ExactMass","Ingredient_smiles")]
      colnames(hmdbdata2) <- c("Compound ID","Metabolites","Formula","Exact_Mass","smiles")
      
      hmdbdata3 <- rbind(hmdbdata1,hmdbdata2)
      
      if(outdatabse == "tcm"){
        
      }else if(outdatabse == "tcme"){
        # tcme <- getmysqldata(dbname = "cosa",
        #                      table = "tcme")
        # hmdbdata3 <- rbind(hmdbdata3,tcme)
      }
      
      hmdbdata3$Formula <- gsub(pattern = "\\+",replacement = "",x = hmdbdata3$Formula)
      hmdbdata3$Formula <- gsub(pattern = "\\-",replacement = "",x = hmdbdata3$Formula)
      hmdbdata3 <- hmdbdata3[!is.na(hmdbdata3$Formula),]
      hmdbdata3 <- hmdbdata3[!is.na(hmdbdata3$Metabolites),]
      hmdbdata3 <- hmdbdata3[!is.na(hmdbdata3$`Compound ID`),]
      hmdbdata3 <- hmdbdata3[!is.na(hmdbdata3$Exact_Mass),]
      hmdbdata3 <- hmdbdata3[!is.na(hmdbdata3$smiles),]
      hmdbdata3 <- hmdbdata3[!grepl(pattern = "or",x = hmdbdata3$Formula),]
      hmdbdata3 <- hmdbdata3[!grepl(pattern = "n",x = hmdbdata3$Formula),]
      hmdbdata3 <- hmdbdata3[grepl(pattern = "^C",x = hmdbdata3$Formula),]
      hmdbdata3 <- hmdbdata3[!duplicated(hmdbdata3$`Compound ID`),]
      hmdbdata3$Exact_Mass <- as.numeric(hmdbdata3$Exact_Mass)
      
      hmdbdata4 <- getmysqldata(dbname = "cosa",
                                table = "compound_structure",
                                tablelist = "cid,inchikey as 'Compound ID'",
                                wherename = "inchikey",
                                wheredata = hmdbdata3$`Compound ID`)
      hmdbdata3 <- merge(hmdbdata3,hmdbdata4,by = "Compound ID",all = T)
      hmdbdata5 <- getmysqldata(dbname = "cosa",
                                table = "compound_classification",
                                tablelist = "cid,superclass as 'Super Class',class as Class,subclass as 'Sub Class'",
                                wherename = "cid",
                                wheredata = hmdbdata4$cid)
      hmdbdata3 <- merge(hmdbdata3,hmdbdata5,by = "cid",all = T)
      hmdbdata6 <- getmysqldata(dbname = "cosa",
                                table = "cid_keggid",
                                tablelist = "cid,keggid as KEGG",
                                wherename = "cid",
                                wheredata = hmdbdata4$cid)
      hmdbdata3 <- merge(hmdbdata3,hmdbdata6,by = "cid",all = T)
      
      for (i in 1:dim(hmdbdata3)[1]) {
        
        SMILES <- hmdbdata3[i,"smiles"]
        
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
      
      adduct <- getmysqldata(dbname = "meta",
                             table = "spatialadduct",
                             wherename = "usetype",
                             wheredata = 1)
      adduct <- adduct[adduct$Type != "isotope",]
      
      for ( i in 1:dim(adduct)[1]) {
        hmdbdata6 <- hmdbdata4
        hmdbdata6$charge <- adduct[i,"charge"]
        hmdbdata6$mode <- adduct[i,"Type"]
        hmdbdata6$Adduct <- adduct[i,"Name"]
        hmdbdata6[,"mz"] <- (hmdbdata6$Exact_Mass*adduct[i,"nM"])/adduct[i,"charge"]+adduct[i,"Mass"]
        hmdbdata5 <- rbind(hmdbdata5,hmdbdata6)
      }
      hmdbdata <- hmdbdata5
      hmdbdata <- hmdbdata[,c("Compound ID","Metabolites","Super Class","Class","Sub Class","KEGG",
                              "Formula","Adduct","mz","mode")]
      
    }else{
      hmdbdata <- readdata(outdatabse)
    }
  }else{
    hmdbdata <- getmysqldata(table = "hmdbmzforsptial2",
                             wherename = "Adduct",
                             wheredata = c("M-","M-2H","M+2H","M+NH4"),
                             wherenotin = T)
    # delcpdid <- getLCMSUntargetedDelCpdsforhmdbid()
    # hmdbdata <- hmdbdata[!(hmdbdata$`Compound ID` %in% delcpdid),]
    hmdbdata <- hmdbdata[!grepl(pattern = "K",x = hmdbdata$Formula),]
    hmdbdata <- hmdbdata[!grepl(pattern = "Na",x = hmdbdata$Formula),]
    hmdbdata <- hmdbdata[!grepl(pattern = "Li",x = hmdbdata$Formula),]
    hmdbdata <- hmdbdata[!grepl(pattern = "Ca",x = hmdbdata$Formula),]
    hmdbdata <- hmdbdata[!grepl(pattern = "Cu",x = hmdbdata$Formula),]
    hmdbdata <- hmdbdata[!grepl(pattern = "Mg",x = hmdbdata$Formula),]
    hmdbdata <- hmdbdata[!grepl(pattern = "Co",x = hmdbdata$Formula),]
    hmdbdata <- hmdbdata[!grepl(pattern = "Cr",x = hmdbdata$Formula),]
    hmdbdata <- hmdbdata[!grepl(pattern = "Fe",x = hmdbdata$Formula),]
    hmdbdata <- hmdbdata[grepl(pattern = "^C",x = hmdbdata$Formula),]
  }
  
  hmdbdata <- hmdbdata[!duplicated(hmdbdata[,c("Metabolites","Adduct")]),]
  hmdbdata2 <- hmdbdata[!duplicated(hmdbdata[,c("Formula","Adduct","mode")]),]
  hmdbdata2 <- hmdbdata2[,c("Formula","Adduct","mz","mode")]
  hmdbdata2 <- hmdbdata2[!is.na(hmdbdata2$Formula),]
  hmdbdata2 <- hmdbdata2[hmdbdata2$mz > 70 & hmdbdata2$mz < 2000,]
  hmdbdata_mode <- hmdbdata2[hmdbdata2$mode == mode,]
  
  # 加和离子与同位素形式
  rawadducts <- getmysqldata(table = "spatialadduct",wherename = "usetype",wheredata = 1)
  isotopedata <- rawadducts[rawadducts$Type == "isotope",]
  isotopedata[,"formula"] <- isotopedata[,"Name"]
  isotopeadducts <- isotopedata[,c("Name","formula","Mass")]
  colnames(isotopeadducts) <- c("name","formula","mass")
  adductsdata <- rawadducts[rawadducts$Type == mode,]
  
  # mass2adduct的相似性计算
  # mz <- readrds(filename = paste0(savepeakpath,"referencemz-raw-", mode, ".rds"))
  # d.diff <- massdiff(mz)
  isotopeannot_cor <- NULL
  adductsannot_cor <- NULL
  mse <- NULL
  for ( i in 1:length(filename)) {
    mse_sample <- readimzml(filename = filename[i],
                            mass.range = mass.range,
                            resolution = resolution,
                            units = units,
                            area = T,
                            areards = areards[i])
    mz <- mz(mse_sample)
    d.diff <- massdiff(mz)
    # mse_sample <- mse_sample %>%
    #   mzAlign(ref = mz,tolerance = 0.05,units ="mz") %>%
    #   peakBin(ref = mz,type=c("height"),tolerance = tolerance, units = units) %>%
    #   process()
    
    # mse_sample <- mse_sample %>%
    #   peakBin(ref = mz,type=c("height"),tolerance = tolerance, units = units) %>%
    #   process()
    
    # spectradata <- iData(mse_sample, "intensity")
    # spectradata <- matrix(spectradata,nrow = dim(spectradata)[1], ncol = dim(spectradata)[2])
    # base::print(dim(spectradata))
    # spectradata2 <- t(matter::apply(spectradata,MARGIN = 1,smooth.image.gaussian,window=3))
    # spectradata2[is.na(spectradata2)] <- 0
    spectradata2 <- spatialsmooth(mse_sample)
    # base::print(dim(spectradata2))
    # base::print((.packages()))
    spectra(mse_sample) <- spectradata2
    
    if(is.null(mse)){
      mse <- mse_sample
    }else{
      mse <- BiocGenerics::cbind(mse, mse_sample)
    }
    
    d <- cardinal2msimat2(mse_sample)
    
    # 同位素
    isotopeannot <- adductMatchnew2(d.diff,add = isotopeadducts,ppm = NULL,mDa = 2)
    if("C2H4O" %in% isotopeadducts$name){
      isotopeannot2 <- adductMatchnew2(d.diff,add = isotopeadducts[isotopeadducts$name == "C2H4O",],ppm = NULL,mDa = 5)
      isotopeannot <- rbind(isotopeannot,isotopeannot2)
      isotopeannot <- isotopeannot[!duplicated(isotopeannot[,c("A","B","matches")]),]
    }
    intensitydataA <- data.frame(mz = mz,
                                 Aintensity = featureApply(mse_sample,mean))
    intensitydataB <- data.frame(mz = mz,
                                 Bintensity = featureApply(mse_sample,mean))
    isotopeannot <- merge(x = isotopeannot,y = intensitydataB,by.x = "B",by.y = "mz")
    isotopeannot <- merge(x = isotopeannot,y = intensitydataA,by.x = "A",by.y = "mz")
    isotopeannot[,"intensityratio"] <- isotopeannot[,"Bintensity"]/isotopeannot[,"Aintensity"]
    isotopeannot <- isotopeannot[!is.na(isotopeannot$intensityratio),]
    isotopeannot <- isotopeannot[(isotopeannot$intensityratio < 0.6) | (isotopeannot$matches == "C2H4O"),]
    j <- 1
    while (j <= dim(isotopeannot)[1]) {
      if(!("C2H4O" %in% strsplit(isotopeannot$matches[j],split = "\\+")[[1]])){
        isotopeannot2 <- isotopeannot[isotopeannot$A %in% isotopeannot$B[j],]
        isotopeannot2 <- isotopeannot2[isotopeannot2$matches != "C2H4O",]
        isotopeannot2 <- isotopeannot2[lapply(strsplit(isotopeannot2$matches,split = "\\+"),length) < (5-length(strsplit(isotopeannot$matches[j],"\\+")[[1]])),]
        if(dim(isotopeannot2)[1] > 0){
          isotopeannot2$matches <- paste0(isotopeannot$matches[j],"+",isotopeannot2$matches)
          isotopeannot2$mass <- isotopeannot$mass[j]+isotopeannot2$mass
          isotopeannot2$A <- isotopeannot$A[j]
          isotopeannot2$diff <- isotopeannot2$diff+isotopeannot$diff[j]
          isotopeannot <- rbind(isotopeannot,isotopeannot2)
        }
        isotopeannot <- isotopeannot[isotopeannot$mass < 1,]
      }
      j <- j+1
    }
    isotopeannot <- isotopeannot[,c("A","B","diff","delta","matches","mass")]
    isotopeannot_cor_single <- corrPairsMSI2(d,isotopeannot)
    isotopeannot_cor <- rbind(isotopeannot_cor,isotopeannot_cor_single)
    
    # 加和离子
    adductsannot <- adductMatchnew(d.diff,adductsdata = adductsdata,ppm = ppm*2)
    adductsannot_cor_single <- corrPairsMSI2(d,adductsannot)
    adductsannot_cor <- rbind(adductsannot_cor,adductsannot_cor_single)
    
    gc(verbose = T)
  }
  
  # 同位素结果
  isotopeannot_cor_1 <- isotopeannot_cor
  isotopeannot_cor_1 <- isotopeannot_cor_1[!is.na(isotopeannot_cor_1$Significance),]
  isotopeannot_cor_1[isotopeannot_cor_1$Estimate < 0,"Significance"] <- 0
  isotopeannot_cor_1[isotopeannot_cor_1$P.value > 0.05,"Significance"] <- 0
  isotopeannot_cor_1 <- isotopeannot_cor_1[(isotopeannot_cor_1$Significance == 1) | (isotopeannot_cor_1$matches == "C2H4O"),]
  isotopeannot_cor_1 <- isotopeannot_cor_1[order(isotopeannot_cor_1$Estimate,decreasing = T),]
  isotopeannot_cor_1 <- isotopeannot_cor_1[order(isotopeannot_cor_1$P.value),]
  isotopeannot_cor_1 <- isotopeannot_cor_1[!duplicated(isotopeannot_cor_1[,c("A","B","matches")]),]
  intensitydataA <- data.frame(mz = mz,
                               Aintensity = featureApply(mse,mean))
  intensitydataB <- data.frame(mz = mz,
                               Bintensity = featureApply(mse,mean))
  isotopeannot_cor_1 <- merge(x = isotopeannot_cor_1,y = intensitydataB,by.x = "B",by.y = "mz")
  isotopeannot_cor_1 <- merge(x = isotopeannot_cor_1,y = intensitydataA,by.x = "A",by.y = "mz")
  isotopeannot_cor_1[,"intensityratio"] <- isotopeannot_cor_1[,"Bintensity"]/isotopeannot_cor_1[,"Aintensity"]
  isotopeannot_cor_1 <- isotopeannot_cor_1[(isotopeannot_cor_1$intensityratio <= 0.6) | (isotopeannot_cor_1$matches == "C2H4O"),]
  
  # isotopedata_abundance <- isotopedata[,c("Name","Abundance")]
  # isotopedata_abundance[,"Abundance"] <- isotopedata_abundance[,"Abundance"]/100
  # isotopeannot_cor_1 <- merge(isotopeannot_cor_1,isotopedata_abundance,by.x = "matches",by.y = "Name")
  
  # 保存同位素结果
  savexlsx1(data = isotopeannot_cor_1,
            filename = paste0(saveannopath,"isotope.xlsx"),
            sheet = mode)
  
  # 加合离子结果
  adductsannot_cor_1 <- adductsannot_cor
  adductsannot_cor_1[is.na(adductsannot_cor_1$Significance),"Estimate"] <- -1
  adductsannot_cor_1[is.na(adductsannot_cor_1$Significance),"P.value"] <- 1
  adductsannot_cor_1[is.na(adductsannot_cor_1$Significance),"Significance"] <- 0
  adductsannot_cor_1 <- adductsannot_cor_1[order(adductsannot_cor_1$Estimate,decreasing = T),]
  adductsannot_cor_1 <- adductsannot_cor_1[order(adductsannot_cor_1$P.value),]
  adductsannot_cor_1 <- adductsannot_cor_1[!duplicated(adductsannot_cor_1[,c("A","B","matches")]),]
  adductsannot_cor_1 <- tidyr::separate(data = adductsannot_cor_1,col = "matches",into = c("matchesA","matchesB"),sep = ";")
  adductsannot_cor_A <- adductsannot_cor_1[,c("A","matchesA","Estimate","P.value","Significance","B","matchesB")]
  colnames(adductsannot_cor_A) <- c("mz","matches","Estimate","P.value","Significance","tomz","tomatches")
  adductsannot_cor_B <- adductsannot_cor_1[,c("B","matchesB","Estimate","P.value","Significance","A","matchesA")]
  colnames(adductsannot_cor_B) <- c("mz","matches","Estimate","P.value","Significance","tomz","tomatches")
  adductsannot_cor_new <- rbind(adductsannot_cor_A,adductsannot_cor_B)
  
  ### 聚合物识别
  if(mode == "neg"){
    polymermz <- NULL
  }else{
    isotopelistdata <- isotopeannot_cor_1[isotopeannot_cor_1$matches == "C2H4O",]
    isotopelistdata <- isotopelistdata[isotopelistdata$A > polymerminmz,]
    polymeradduct <- adductsannot_cor_new[adductsannot_cor_new$Estimate > 0,]
    polymeradduct <- polymeradduct[polymeradduct$mz > polymerminmz,]
    polymeradduct <- polymeradduct[polymeradduct$tomz > polymerminmz,]
    
    if(dim(isotopelistdata)[1] == 0){
      polymermz <- NULL
    }else{
      isotopelistdata <- isotopelistdata[order(isotopelistdata$A),]
      isotopelistdata[,"type"] <- NA
      
      for( i in 1:dim(isotopelistdata)[1]){
        isotopelistdataA <- isotopelistdata[isotopelistdata[,"B"] == isotopelistdata[i,"A"],]
        if(dim(isotopelistdataA)[1] == 0){
          isotopelistdata[i,"type"] <- RandomCode()
        }else{
          isotopelistdata[i,"type"] <- isotopelistdataA$type[1]
        }
      }
      
      isomzdata <- table(isotopelistdata$type)
      isomzdata <- isomzdata[order(isomzdata,decreasing = T)]
      isomzdata <- isomzdata[isomzdata >= polymernum]
      isotopelistdata <-  isotopelistdata[isotopelistdata$type %in% names(isomzdata),]
      polymeradduct <- polymeradduct[polymeradduct$mz %in% c(isotopelistdata$A,isotopelistdata$B),]
      
      isotopelistdata2 <- isotopeannot_cor_1[isotopeannot_cor_1$matches != "C2H4O",]
      isotopelistdata2 <- isotopelistdata2[isotopelistdata2$A %in% c(isotopelistdata$A,isotopelistdata$B,polymeradduct$tomz),]
      polymermz <- c(isotopelistdata$A,isotopelistdata$B,polymeradduct$tomz,isotopelistdata2$B)
      polymermz <- unique(polymermz)
      polymermz <- polymermz[order(polymermz)]
      
      print(paste0("聚合物数量:",length(polymermz)))
      
      saverds(data = polymermz,
              filename = paste0(saveannopath,"polymermz-", mode, ".rds"))
      
    }
  }
  
  mzlist <- mz(mse)
  mzlist <- mzlist[!(mzlist %in% polymermz)]
  adductsannot_cor_new2 <- adductsannot_cor_new[!(adductsannot_cor_new$mz %in% polymermz),]
  adductsannot_cor_new2 <- adductsannot_cor_new2[!(adductsannot_cor_new2$tomz %in% polymermz),]
  mzlist <<- mzlist
  if(waters){
    print("~使用waters定性模式")
    hmdbdata_mode_mz <- getmatchformula(mzdata = mzlist,formuladata = hmdbdata_mode,mode = mode,ppm = ppm*2)
    hmdbdata_mode_mz <- hmdbdata_mode_mz[(!(hmdbdata_mode_mz$ppm > ppm)) | hmdbdata_mode_mz$mz < 400,]
    hmdbdata_mode_mz <<- hmdbdata_mode_mz
  }else{
    hmdbdata_mode_mz <- getmatchformula(mzdata = mzlist,formuladata = hmdbdata_mode,mode = mode,ppm = ppm)
  }
  adductsannot_cor_new2 <- merge(adductsannot_cor_new2,hmdbdata_mode_mz,by.x = c("mz","matches"),by.y=c("mz","Adduct"))
  colnames(hmdbdata_mode_mz)[2] <- "matches"
  hmdbdata_mode_mz[,c("Estimate","P.value","Significance","tomz","tomatches")] <- NA
  adductsannot_cor_new2 <- rbind(adductsannot_cor_new2,hmdbdata_mode_mz)
  adductsannot_cor_new2[is.na(adductsannot_cor_new2$Estimate),"Estimate"] <- 0
  if(mode == "neg"){
    # adductsannot_cor_new2[adductsannot_cor_new2$matches == "M-","Estimate"] <- 0.1
  }else{
    adductsannot_cor_new2[adductsannot_cor_new2$matches == "M+","Estimate"] <- 0.1
  }
  
  adductsannot_cor_new3 <- adductsannot_cor_new2
  # adductsannot_cor_new3 <- adductsannot_cor_new2 %>% 
  #   group_by(mz,matches,Formula) %>% 
  #   mutate(adductscore = mean(Estimate)) %>%
  #   select(!c(tomz, tomatches)) %>%
  #   distinct(mz,matches,.keep_all = T)
  # adductsannot_cor_new3 <- as.data.frame(adductsannot_cor_new3)
  adductsdata2 <- adductsdata[,c("Name","Element","nM","importance")]
  colnames(adductsdata2) <- c("matches","adductElement","nM","adductimportance")
  adductsannot_cor_new3 <- merge(adductsannot_cor_new3,adductsdata2,by = "matches",all.x = T)
  
  for ( i in 1:dim(adductsannot_cor_new3)[1]) {
    if(is.na(adductsannot_cor_new3[i,"adductElement"])){
      adductsannot_cor_new3[i,"Formula-adduct"] <- adductsannot_cor_new3[i,"Formula"]
    }else{
      adductsannot_cor_new3[i,"Formula-adduct"] <- getnewformula(formula = adductsannot_cor_new3[i,"Formula"],
                                                                 adductElement = adductsannot_cor_new3[i,"adductElement"],
                                                                 nM = adductsannot_cor_new3[i,"nM"])
    }
  }
  
  adductsannot_cor_new3 <- adductsannot_cor_new3[!is.na(adductsannot_cor_new3[,"Formula-adduct"]),]
  isotopeannot_cor_2 <- isotopeannot_cor_1[isotopeannot_cor_1$matches != "C2H4O",]
  isotopeannot_cor_2 <- isotopeannot_cor_2[!(isotopeannot_cor_2$A %in% polymermz),]
  isotopeannot_cor_2 <- isotopeannot_cor_2[!(isotopeannot_cor_2$B %in% polymermz),]
  colnames(isotopeannot_cor_2) <- c("mz","isomz","isomatch","isodiff","isoscore","isop.vale",
                                    "isoSignificance","isointensity","mzintensity","intensityratio")
  adductsannot_cor_new4 <- merge(adductsannot_cor_new3,isotopeannot_cor_2,by = "mz")
  adductsannot_cor_new5 <- data.frame(adductsannot_cor_new3,
                                      "isomz" = NA,"isomatch" = NA,"isodiff" = NA,"isoscore" = 0,"isop.vale" = NA,
                                      "isoSignificance" = NA,"isointensity" = NA,"mzintensity" = NA,"intensityratio" = NA,
                                      check.names = F)
  adductsannot_cor_new4 <- rbind(adductsannot_cor_new5,adductsannot_cor_new4)
  for ( i in 1:dim(adductsannot_cor_new4)[1]) {
    if(is.na(adductsannot_cor_new4$isodiff[i])){
      adductsannot_cor_new4[i,"select"] <- T
      next
    }
    adductsannot_cor_new4[i,"select"] <-  all(unlist(stringr::str_extract_all(string = adductsannot_cor_new4$isodiff[i],pattern = "[A-Z][a-z]*")[[1]]) %in% 
                                                unlist(stringr::str_extract_all(string = adductsannot_cor_new4$`Formula-adduct`[i],pattern = "[A-Z][a-z]*")[[1]]))
  }
  adductsannot_cor_new4 <- adductsannot_cor_new4[order(adductsannot_cor_new4$mz),]
  adductsannot_cor_new4 <- adductsannot_cor_new4[order(adductsannot_cor_new4$isoscore,decreasing = T),]
  adductsannot_cor_new4 <- adductsannot_cor_new4[order(adductsannot_cor_new4$select,decreasing = T),]
  adductsannot_cor_new4 <- adductsannot_cor_new4[!duplicated(adductsannot_cor_new4[,c("matches","mz","Formula","isomz","tomz")]),]
  adductsannot_cor_new4[!adductsannot_cor_new4$select,"isoscore"] <- 0
  adductsannot_cor_new4[is.na(adductsannot_cor_new4$adductimportance),"adductimportance"] <- 0
  # adductsannot_cor_new4 <- adductsannot_cor_new4 %>% 
  #   group_by(mz,matches,Formula) %>% 
  #   mutate(isoscoremean = mean(isoscore))
  # adductsannot_cor_new4 <- as.data.frame(adductsannot_cor_new4[order(adductsannot_cor_new4$matches),])
  # adductsannot_cor_new4[is.na(adductsannot_cor_new4$adductimportance),"adductimportance"] <- 0
  # adductsannot_cor_new4[,"sumscore"] <- adductsannot_cor_new4$adductscore/10*7+
  #   (ppm-adductsannot_cor_new4$ppm)/ppm/20+
  #   adductsannot_cor_new4$isoscoremean/10*2+
  #   (max(adductsannot_cor_new4$adductimportance)-adductsannot_cor_new4$adductimportance)/max(adductsannot_cor_new4$adductimportance)/20
  # adductsannot_cor_new4 <- adductsannot_cor_new4[order(adductsannot_cor_new4$sumscore,decreasing = T),]
  adductsannot_cor_new5 <- adductsannot_cor_new4
  adductsannot_cor_new6 <- adductsannot_cor_new4[0,]
  
  meancal <- function(x){
    x <- x[!duplicated(x)]
    x <- x[!is.na(x)]
    x <- x[x!=0]
    if(length(x) == 0){
      return(0)
    }else{
      return(mean(x))
    }
  }
  meancal2 <- function(x){
    x <- x[!duplicated(x)]
    x <- x[!is.na(x)]
    if(length(x) == 0){
      return(0)
    }else{
      return(mean(x))
    }
  }
  
  while (dim(adductsannot_cor_new5)[1] != 0) {
    
    adductsannot_cor_new5 <- adductsannot_cor_new5 %>%
      group_by(mz,matches,Formula,isodiff) %>%
      mutate(adductscore = meancal(Estimate)) %>%
      group_by(mz,matches,Formula,isodiff) %>%
      mutate(isoscoremean = meancal2(isoscore))
    adductsannot_cor_new5 <- as.data.frame(adductsannot_cor_new5)
    adductsannot_cor_new5[,"sumscore"] <- adductsannot_cor_new5$adductscore/10*7+
      (ppm-adductsannot_cor_new5$ppm)/ppm/20+
      adductsannot_cor_new5$isoscoremean/10*2+
      (max(adductsannot_cor_new5$adductimportance)-adductsannot_cor_new5$adductimportance)/max(adductsannot_cor_new5$adductimportance)/20
    adductsannot_cor_new5 <- adductsannot_cor_new5[order(adductsannot_cor_new5$sumscore,decreasing = T),]
    
    adductsannot_cor_new5_1 <- adductsannot_cor_new5[adductsannot_cor_new5$Formula == adductsannot_cor_new5$Formula[1],]
    adductsannot_cor_new5_1 <- adductsannot_cor_new5_1[adductsannot_cor_new5_1$select,]
    adductsannot_cor_new5_1 <- adductsannot_cor_new5_1[adductsannot_cor_new5_1$adductscore > mean(unique(adductsannot_cor_new5_1$adductscore))-0.1,]
    if(all(adductsannot_cor_new5_1$adductscore <= 0)){
      adductsannot_cor_new5_1 <- adductsannot_cor_new5_1[adductsannot_cor_new5_1$ppm == min(adductsannot_cor_new5_1$ppm),]
    }else{
      adductsannot_cor_new5_1 <- adductsannot_cor_new5_1[adductsannot_cor_new5_1$adductscore > 0,]
    }
    imagemz <- unique(c(adductsannot_cor_new5_1$mz,adductsannot_cor_new5_1$isomz))
    imagemz <- imagemz[!is.na(imagemz)]
    # image(mse,mz = imagemz, normalize.image="linear")
    adductsannot_cor_new6 <- rbind(adductsannot_cor_new6,adductsannot_cor_new5_1)
    adductsannot_cor_new5 <- adductsannot_cor_new5[!(adductsannot_cor_new5$mz %in% imagemz),]
    adductsannot_cor_new5 <- adductsannot_cor_new5[!(adductsannot_cor_new5$isomz %in% imagemz),]
    adductsannot_cor_new5 <- adductsannot_cor_new5[!(adductsannot_cor_new5$tomz %in% imagemz),]
    adductsannot_cor_new5 <- adductsannot_cor_new5[!(adductsannot_cor_new5$Formula %in% adductsannot_cor_new5_1$Formula),]
  }
  
  savexlsx1(data = adductsannot_cor_new4,
            filename = paste0(saveannopath,"adducts-mid.xlsx"),
            sheet = mode)
  
  # mz	adducts	intensity	value	type	score	isotope	formula	sumscore	ppmscore	ppm	polymer
  mzmetadata_1 <-  adductsannot_cor_new6[,c("mz","matches","Formula","ppm","sumscore")]
  colnames(mzmetadata_1) <- c("mz","adducts","formula","ppm","score")
  mzmetadata_1$isotope <- NA
  mzmetadata_2 <-  adductsannot_cor_new6[,c("isomz","matches","Formula","ppm","sumscore","isodiff")]
  colnames(mzmetadata_2) <- c("mz","adducts","formula","ppm","score","isotope")
  mzmetadata <- rbind(mzmetadata_1,mzmetadata_2)
  mzmetadata <- mzmetadata[!is.na(mzmetadata$mz),]
  mzmetadata <- mzmetadata[!duplicated(mzmetadata$mz),]
  
  rawmzdata <- data.frame(mz = mz,intensity = featureApply(mse,mean))
  mzmetadata2 <- merge(rawmzdata,mzmetadata,by = "mz",all.x = T,sort = F)
  mzmetadata2 <- mzmetadata2[order(mzmetadata2$mz),]
  mzmetadata2[mzmetadata2$mz %in% polymermz,"polymer"] <- "C2H4O"
  
  saverds(data =  mzmetadata2,
          filename = paste0(saveannopath,"adducts-", mode, ".rds"))
  savexlsx1(data = mzmetadata2,
            filename = paste0(saveannopath,"adducts.xlsx"),
            sheet = mode)
  
  hmdbdata <- hmdbdata[hmdbdata$mode == mode,]
  hmdbdata3 <- hmdbdata[!duplicated(hmdbdata$`Compound ID`),]
  hmdbdata3 <- hmdbdata3[order(hmdbdata3$`Compound ID`),]
  hmdbdata3 <- hmdbdata3[order(hmdbdata3$KEGG),]
  hmdbdata3 <- hmdbdata3[!duplicated(hmdbdata3[,c("Metabolites","Formula")]),]
  
  mzmetadata3 <- mzmetadata2[!is.na(mzmetadata2$adducts),]
  mzmetadata3 <- mzmetadata3[!duplicated(mzmetadata3$formula),]
  
  qualitativedata <- NULL
  
  for ( i in 1:dim(mzmetadata3)[1]) {
    hmdbdata4 <- hmdbdata3[hmdbdata3$Formula == mzmetadata3$formula[i],]
    hmdbdata4 <- hmdbdata4[,c("Compound ID","Metabolites","Super Class","Class","Sub Class","KEGG","Formula","mode")]
    hmdbdata4[,"level"] <- NA
    hmdbdata4[!is.na(hmdbdata4$KEGG),"level"] <- "*"
    qualitativedata <- rbind(qualitativedata,hmdbdata4)
  }
  
  species_tissuedata <- getspecies_tissuedata(species_tissue = species_tissue,
                                              copypath = saveannopath)
  lcqualitativefrom <- species_tissuedata$lcqualitativefrom
  manqualitativefrom <- species_tissuedata$manqualitativefrom
  species_tissuelocation <- species_tissuedata$hmdbloaction
  
  if (file.exists(lcqualitativefrom)) {
    mandata1 <- readdata(filename = lcqualitativefrom,sheet = "pos")
    mandata2 <- readdata(filename = lcqualitativefrom,sheet = "neg")
    mandata <- rbind(mandata1,mandata2)
    mandata <- mandata[!duplicated(mandata$id),]
    
    qualitativedata[qualitativedata$`Compound ID` %in% mandata$id,"level"] <- "**"
  }else{
    print(paste0(lcqualitativefrom, "文件不存在"))
  }
  
  
  if (file.exists(manqualitativefrom)) {
    mandata1 <- readdata(filename = manqualitativefrom,sheet = "pos")
    mandata2 <- readdata(filename = manqualitativefrom,sheet = "neg")
    mandata <- rbind(mandata1,mandata2)
    mandata <- mandata[!duplicated(mandata$`Compound ID`),]
    
    qualitativedata[qualitativedata$`Compound ID` %in% mandata$`Compound ID`,"level"] <- "***"
  }else{
    print(paste0(manqualitativefrom, "文件不存在"))
  }
  
  qualitativedata2 <- qualitativedata %>% group_by(Formula) %>% mutate(num = length(Formula))
  
  savexlsx1(data = qualitativedata2,
            filename = paste0(saveannopath,"anno-all.xlsx"),
            sheet = mode)
  
  formuladata <- qualitativedata$Formula
  formuladata <- formuladata[!duplicated(formuladata)]
  
  qualitativedata3 <- NULL
  for ( i in 1:length(formuladata)) {
    qualitativedata2 <- qualitativedata[qualitativedata$Formula == formuladata[i],]
    if(any(qualitativedata2$level %in% "***")){
      qualitativedata2 <- qualitativedata2[qualitativedata2$level %in% "***",]
    }else if(any(qualitativedata2$level %in% "**")){
      qualitativedata2 <- qualitativedata2[qualitativedata2$level %in% "**",]
    }else if(any(qualitativedata2$level %in% "*")){
      qualitativedata2 <- qualitativedata2[qualitativedata2$level %in% "*",]
    }
    
    if(dim(qualitativedata2)[1] > 1){
      if(any(!is.na(qualitativedata2$KEGG))){
        qualitativedata2 <- qualitativedata2[!is.na(qualitativedata2$KEGG),]
      }
      
      if(dim(qualitativedata2)[1] > 3){
        if(any(qualitativedata2$`Compound ID` %in% species_tissuelocation)){
          qualitativedata2 <- qualitativedata2[qualitativedata2$`Compound ID` %in% species_tissuelocation,]
        }
      }
      
      if(dim(qualitativedata2)[1] > 3){
        qualitativedata2 <- qualitativedata2[order(qualitativedata2$`Compound ID`),]
        qualitativedata2 <- qualitativedata2[1:3,]
      }
    }
    
    qualitativedata3 <- rbind(qualitativedata3,qualitativedata2)
  }
  
  qualitativedata3 <- qualitativedata3 %>% group_by(Formula) %>% mutate(num = length(Formula))
  
  savexlsx1(data = qualitativedata3,
            filename = paste0(saveannopath,"anno.xlsx"),
            sheet = mode)
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-i","--imzmlpath",default = "./sample/imzml-pre/", help = "imzml原始文件路径,默认./sample/imzml-pre/")
  parser$add_argument("-a","--areapath",default = "./sample/select/", help = "选区文件路径,默认./sample/select/")
  parser$add_argument("-sp","--savepeakpath",default = "./sample/peak/", help = "mz数据保存路径,默认./sample/peak/")
  parser$add_argument("-m","--massrange",default = NULL,type= "double",help = "mz范围,默认NULL",
                      nargs = 2,dest = "mass.range")
  parser$add_argument("-re","--resolution",default = 5, type= "double",help = "分辨率,默认5")
  parser$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  parser$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  
  parser$add_argument("-sa","--sample",default = NULL,help = "选择对齐样本")
  parser$add_argument("-ds","--delsample",default = c("bg_data","qc_data"),help = "删除对齐样本",nargs = "+")
  
  # 定性参数
  parser$add_argument("-tr","--tolerance",default = 5, type= "double",help = "峰对齐参数,默认5")
  parser$add_argument("-p","--ppm",default = 5, type= "double",help = "定性ppm参数")
  parser$add_argument("-pn","--polymernum",default = 3, type= "integer",help = "胶聚合物最小数量")
  parser$add_argument("-pm","--polymerminmz",default = 300, type= "double",help = "聚合物最小mz")
  parser$add_argument("-st","--species_tissue",default = "小鼠",help = "物种_组织",required = T)
  parser$add_argument("-sap","--saveannopath",default = "./sample/qualitative/", help = "定性保存路径,默认./sample/qualitative/")
  parser$add_argument("-od","--outdatabse",default = NULL,help = "外部数据库，如是中药数据库输入tcm")
  parser$add_argument("-w","--waters",default = F,help = "是否使用waters仪器",action ='store_true')
  
  parser$add_argument("-mn","--maxnum",default = 10, type= "integer",help = "最大样本数")
  parser$add_argument("-se","--seed",default = 1111, type= "integer",help = "样本选择随机种子")
  
  args <- parser$parse_args()
  writeinfo()
  
  createdir(filename = args$saveannopath,linkdir = T)
  result <- do.call(what = imzmlpreQualitative,args = args)
  
  writeinfo(endtime = T)
}
