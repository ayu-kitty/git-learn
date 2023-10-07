#!/opt/conda/bin/Rscript

#' datafetch_untargetgcms_obj
#'
#' GC非靶数据获取
#'
#' @param data obj
#' @param ... 见`datafetch_untargetgcms`
#'
#' @export
datafetch_untargetgcms_obj_2 <- function(data,mode="fame",
                                         ...){
  if(is.na(grep("细胞",data$info$basic$样本类型)[1])){
    data$data$rawdata$data <- datafetch_untargetgcms_2(filename = data$info$datafrom$GC原始文件,mode = mode,
                                                     ...)
  }else{
    data$data$rawdata$data <- datafetch_untargetgcms_2(filename = data$info$datafrom$GC原始文件,
                                                       mode = "famearea",...)
  }
  
  
  return(data)
}

#' datafetch_untargetgcms
#' 
#' 获取gc定性定量文件2023.7版本
#' a)峰面积归一化算法:由 样本峰面积/总峰面积*10000 修改为 样本峰面积/总峰面积*所有样本总峰面积均值
#' b)内标分段归一化算法:由 分段样本峰面积/分段内标峰面积 修改为分段样本峰面积/分段内标峰面积* 分段所有样本内标峰面积均值
#' c)famearea细胞样本峰面积+内标归一化(先进行内标分段归一化: 分段样本峰面积/分段内标峰面积* 分段所有样本内标峰面积均值;再进行峰面积归一化: 样本峰面积/总峰面积*所有样本总峰面积均值(内标归一化后))
#'
#' @export 
datafetch_untargetgcms_2 <- function(filename = NULL,
                                     samplename = NULL,
                                     grouping = NULL,
                                     mode = "fame", # area峰面积,fame内标分段归一化,single单一内标归一化,famearea细胞样本峰面积+内标归一化
                                     projectclass = "gc", # spme(顶空) or gc(非靶gc) or wax(蜡质)
                                     famersd = 0.1, # 内标rsd筛选标准0.15
                                     score = 70,
                                     savefile = "GC-预处理矩阵.xlsx",
                                     fameaverage = F,
                                     RI = T,
                                     peaklimit = 5000,
                                     singlefame = "Decanoic Acid, Ethyl Ester", # 单内标名字Decanoic acid, ethyl ester葵酸乙酯
                                     deallist = switch (projectclass,
                                                        "gc" = c("Fame C", "Unknown", "Tbs Compound", "Z [^a-z]", 
                                                                 "3,4-Dichlorophenylalanine","fluo", "chlo", 
                                                                 "brom", "iod", "sulf", "silane","sil"),
                                                        "spme" = c("Fame C", "Unknown", "Tbs Compound", "Z [^a-z]", 
                                                                   "3,4-Dichlorophenylalanine","fluo", "chlo", 
                                                                   "brom", "iod", "sulf", "silane","sil",
                                                                   "tms derivative", "Pyridine", "Propoxur", 
                                                                   "Ethanol", "tbdms derivative"),
                                                        "wax" = c("Fame C", "Unknown", "Tbs Compound", "Z [^a-z]", 
                                                                  "3,4-Dichlorophenylalanine","fluo", "chlo", 
                                                                  "brom", "iod", "sulf", "silane","sil",
                                                                  "pyridine", "propoxur", 
                                                                  "ethanol", "tbdms derivative"))){
  suppressMessages(library("dplyr"))
  
  datainfo1 <- data.frame(filename = filename, 
                          mode = mode, 
                          projectclass = projectclass,
                          famersd = famersd,
                          score = score)
  datainfo2 <- data.frame("deallist" = deallist)
  
  # 处理参数保存
  savexlsx1(filename = savefile, sheet = "运行参数", data = datainfo1)
  savexlsx1(filename = savefile, sheet = "剔除列表", data = datainfo2)
  
  # 数据读入
  spmedata <- readdata(filename = filename, header = F, skip = 0)
  famenorm <- data.frame()
  
  # 数据整理
  colnames(spmedata) <- spmedata[5, ]
  colnames(spmedata)[5] <- "Metabolites"
  spmedata$Metabolites <- stringr::str_to_title(spmedata$Metabolites)
  # 删除空的表头
  dataraw <- spmedata[6:dim(spmedata)[1], 1:max(which(spmedata[2, ] == "Sample" | spmedata[2, ] == "QC"))]
  rownames(dataraw) <- dataraw[, 1]
  # 定性数据提出
  dataID <- dataraw[, c("Alignment ID","Average Rt(min)","Quant mass","Metabolites","INCHIKEY","Total score","Total spectrum similarity")]
  dataID[dataID[,6]=="null",6]<-NA
  dataID[dataID[,7]=="null",7]<-NA
  dataID[, c(1:3,6,7)] <- as.numeric(unlist(dataID[, c(1:3,6,7)]))
  # 定量数据提出
  dataM <- dataraw[, -1:-28]
  # base::print(colnames(dataM))
  dataM[,] <- apply(dataraw[, -1:-28], 2,as.numeric)
  # 定量矩阵名字排序和替换为samplename
  if(!is.null(samplename)){
    if (length(samplename) == length(grouping) & all(grouping %in% colnames(dataM))) {
      dataM <- dataM[, grouping]
      colnames(dataM) <- samplename
    } else {
      print("样本量与样品信息不符合")
      stop("样本量与样品信息不符合")
    } 
  }
  #低于5000的峰归0
  dataM[dataM<= peaklimit] <- 0
  # 原始数据保存
  dataraw <- merge(dataID, dataM, by=0,sort = F)
  rownames(dataraw) <- dataraw[,1]
  dataraw <- dataraw[,-1]
  savexlsx1(filename = savefile, sheet = "原始数据矩阵", data = dataraw)
  # 内标数据提取
  dataFameC <- dataraw[grep("^Fame C", dataraw$Metabolites), ]
  dataFameCM <- dataFameC[, -1:-6]
  
  singlefame_area <- dataraw[grep(singlefame, dataraw$Metabolites), ]
  singlefame_area <- arrange(singlefame_area, desc(`Total score`))
  singlefame_area <- singlefame_area[!duplicated(singlefame_area$Metabolites),]
  singlefameM <- singlefame_area[,-1:-6]
  
  # 内标数据0值检查及处理
  if (mode == "fame") {
    
    if(dim(dataFameCM)[1] == 0){
      stop("Fame内标异常:未找到任何内标")
    }
    
    if (fameaverage) {
      print("使用fameaverage = T,处理内标中")
      for (i in 1:dim(dataFameCM)[1]) {
        dataFameCM[i,dataFameCM[i, ]==0] <- mean(t(dataFameCM[i,dataFameCM[i, ]!=0]))
      }
      #2022.11.25新增传参
      dataFameC[,7:length(dataFameC)] <- dataFameCM
    }
    
    n <- sum(dataFameCM == 0)
    if (n != 0) {
      warning(paste0("Fame内标异常:", n, "个内标为空值,请检查内标!!!检查无误后请使用fameaverage = T,后续将以该内标所有样本均值替代空值"))
      for (i in 1:dim(dataFameCM)[1]) {
        dataFameCM[i,dataFameCM[i, ]==0] <- mean(t(dataFameCM[i,dataFameCM[i, ]!=0]))
      }
      #2022.11.25新增传参
      dataFameC[,7:length(dataFameC)] <- dataFameCM
    } else {
      print("内标无异常")
    }
  } else if (mode == "single") {
    if(dim(singlefameM)[1] == 0){
      stop("single内标异常:未找到任何内标")
    }
    
    if (fameaverage) {
      print("使用fameaverage = T,处理内标中")
      for (i in 1:dim(singlefameM)[1]) {
        singlefameM[i, singlefameM[i,]==0] <- mean(t(singlefameM[i,singlefameM[i, ]!=0]))
      }
    }
    
    n <- sum(singlefameM == 0)
    if (n != 0) {
      stop(paste0("single内标异常:", n, "个内标为空值,请检查内标!!!检查无误后请使用fameaverage = T,后续将以该内标所有样本均值替代空值"))
    } else {
      print("内标无异常")
    }
  }
  
  # 峰面积归一化
  if (mode == "area") {
    
    print("当前正在使用的归一化方法:峰面积归一化")
    # 计算总峰面积
    allarea <- apply(X = dataM, MARGIN = 2, FUN = sum)
    # 归一化后ID和M文件合并并保存
    data1 <- data.frame(dataID, t(t(dataM) / allarea * mean(allarea)), check.names = F)
    savexlsx1(filename = savefile, sheet = "峰面积归一化", data = data1)
    
  } else if (mode == "fame") {
    
    # 内标分段归一化
    print("当前正在使用的归一化方法:内标分段归一化")
    # 计算所有样本内标rsd并保存
    datarsd <- data.frame(
      `mean-ALL` = apply(dataFameC[, 8:length(dataFameC)], MARGIN = 1, FUN = mean),
      `sd-ALL` = apply(dataFameC[, 8:length(dataFameC)], MARGIN = 1, FUN = sd),
      check.names = F
    )
    # 计算qc样本的内标rsd并保存
    # datarsd <- data.frame(`mean-qc` = apply(dataFameC[,grep("QC",colnames(dataFameC),ignore.case = T)], MARGIN = 1, FUN = mean),
    #                      `sd-qc` = apply(dataFameC[,grep("QC",colnames(dataFameC),ignore.case = T)], MARGIN = 1, FUN = sd),
    #                      check.names = F)
    # dataFameC$'rsd-qc' <- datarsd$`sd-qc`/datarsd$`mean-qc`
    dataFameC$"rsd-ALL" <- datarsd$`sd-ALL` / datarsd$`mean-ALL`
    dataFameC <- dataFameC[, c(1:7, length(dataFameC), 8:(length(dataFameC) - 1))]
    savexlsx1(filename = savefile, sheet = "FAME", data = dataFameC)
    
    # 内标分段归一化
    # 内标rsd筛选
    datafamersd <- dataFameC %>% filter(dataFameC$`rsd-ALL` < famersd)
    if(dim(datafamersd)[1] == 0){
      warning("内标rsd均不符合条件")
      datafamersd <- dataFameC %>% filter(dataFameC$`rsd-ALL` < 0.5)
    }
    datafameM <- datafamersd[9:length(datafamersd[1, ])]
    timelist <- c(0, datafamersd$`Average Rt(min)`)
    # 分段归一(part1)
    for (i in 2:(length(timelist))) {
      datasub <- dataM %>% filter(dataraw$`Average Rt(min)` > timelist[i - 1], dataraw$`Average Rt(min)` <= timelist[i])
      famenorm <- rbind(famenorm, t(t(datasub) / as.numeric(datafameM[i - 1, ]) * mean(as.numeric(datafameM[i - 1, ]))))
    }
    # 分段归一(part2)
    datasub <- dataM %>% filter(dataraw$`Average Rt(min)` > timelist[length(timelist)])
    famenorm <- rbind(famenorm, datasub / as.numeric(datafameM[length(datafameM[, 1]), ]) * mean(as.numeric(datafameM[length(datafameM[, 1]), ])))
    # 内标分段归一化结果M与ID合并
    data1 <- data.frame(dataID, famenorm, check.names = F)
    savexlsx1(filename = savefile, sheet = "内标分段归一化", data = data1)
    
  } else if( mode == "single"){
    # 单一内标归一化
    print("当前正在使用的归一化方法:单一内标归一化")
    # 提取单内标
    savexlsx1(filename = savefile,sheet = "singlefame",data = singlefame_area)
    # 归一化后ID和M文件合并并保存
    data1 <- data.frame(dataID, t(t(dataM) / as.numeric(singlefameM) * mean(as.numeric(singlefameM))), check.names = F)
    savexlsx1(filename = savefile, sheet = paste0(singlefame, "归一化"), data = data1)
  }else if (mode == "famearea"){
    print("当前正在使用的归一化方法:内标分段归一化+峰面积归一化")
    # 计算所有样本内标rsd并保存
    datarsd <- data.frame(
      `mean-ALL` = apply(dataFameC[, 8:length(dataFameC)], MARGIN = 1, FUN = mean),
      `sd-ALL` = apply(dataFameC[, 8:length(dataFameC)], MARGIN = 1, FUN = sd),
      check.names = F
    )
    dataFameC$"rsd-ALL" <- datarsd$`sd-ALL` / datarsd$`mean-ALL`
    dataFameC <- dataFameC[, c(1:7, length(dataFameC), 8:(length(dataFameC) - 1))]
    savexlsx1(filename = savefile, sheet = "FAME", data = dataFameC)
    
    # 内标分段归一化
    # 内标rsd筛选
    datafamersd <- dataFameC %>% filter(dataFameC$`rsd-ALL` < famersd)
    if(dim(datafamersd)[1] == 0){
      stop("内标rsd均不符合条件")
    }
    datafameM <- datafamersd[9:length(datafamersd[1, ])]
    timelist <- c(0, datafamersd$`Average Rt(min)`)
    # 分段归一(part1)
    for (i in 2:(length(timelist))) {
      datasub <- dataM %>% filter(dataraw$`Average Rt(min)` > timelist[i - 1], dataraw$`Average Rt(min)` <= timelist[i])
      famenorm <- rbind(famenorm, t(t(datasub) / as.numeric(datafameM[i - 1, ]) * mean(as.numeric(datafameM[i - 1, ]))))
    }
    # 分段归一(part2)
    datasub <- dataM %>% filter(dataraw$`Average Rt(min)` > timelist[length(timelist)])
    famenorm <- rbind(famenorm, datasub / 
                        as.numeric(datafameM[length(datafameM[, 1]), ]) * 
                        mean(as.numeric(datafameM[length(datafameM[, 1]), ])))
    # 计算总峰面积
    allarea2 <- apply(X = famenorm, MARGIN = 2, FUN = sum)
    # 归一化后ID和M文件合并并保存
    data1 <- data.frame(dataID, t(t(famenorm) / allarea2 * mean(allarea2)), check.names = F)
    savexlsx1(filename = savefile, sheet = "内标分段峰面积归一化", data = data1)

  }
  
  # 筛选
  # 1.编号处理
  data1$Metabolites <- gsub(";[0-9]+$", "", data1$Metabolites)
  # 2.删除unknown/z /tbs等代谢物
  data2 <- data1[-(grep(paste0(tolower(deallist), collapse = "|"),tolower(data1$Metabolites))),]
  data2 <- data2[data2$INCHIKEY != "nan",]
  # 3.删除打分低于score代谢物
  data2 <- data2 %>% filter(data2$`Total score` > score)
  # 4.1去除名字重复
  data2 <- arrange(data2, desc(`Total score`))
  data2 <- data2[!duplicated(data2$Metabolites), ]
  # 4.2去除inchikey重复
  data2 <- data2[!duplicated(data2$INCHIKEY), ]
  
  # 数据注释
  data3 <- getmetainfo(data = data2[,"INCHIKEY",drop = F],idlist = "INCHIKEY",
                       needlist = c("cid",
                                    "CAS","HMDB","METLIN","Lipidmaps","KEGG","PubChem","ChEBI",
                                    "smiles","InChIKey","Formula",
                                    "Super Class","Class","Sub Class"))
  
  data3 <- merge(x = data2[,1:6],y = data3,by = "INCHIKEY")
  data3 <- merge(x = data3,y = data2[,c(5,7:dim(data2)[2])],by = "INCHIKEY")[,-1]
  data3 <- data3[order(data3$`Alignment ID`),]
  # data3 <- data3[,-1]
  # data3 <- data3[, c(11:15,2:7,1,8:10,16:dim(data3)[2])]
  savexlsx1(filename = savefile, sheet = "注释结果", data = data3)
  
  # RI指数计算
  if (RI == T) {
    datari <- getmysqldata(dbname = "cosa",
                           table = "gc_ri")
    dataFameC <- merge(x = dataFameC, y = datari, by = "Metabolites")
    dataFameC <- dataFameC[order(dataFameC$`Average Rt(min)`),]
    data3 <- data3[order(data3$`Average Rt(min)`), ]
    i <- 1
    while (i <= dim(data3)[1]) {
      
      j <- which(data3[i, "Average Rt(min)"] > dataFameC$`Average Rt(min)`)
      
      if (length(j) == 0) {
        j <- 1
      } else if (length(j) == dim(dataFameC)[1]) {
        j <- dim(dataFameC)[1] - 1
      } else {
        j <- length(j)
      }
      
      data3[i, "Average RI"] <- dataFameC[j, "RI"] +(dataFameC[j + 1, "RI"] - dataFameC[j, "RI"]) *(data3[i, "Average Rt(min)"] - dataFameC[j, "Average Rt(min)"]) / (dataFameC[j + 1, "Average Rt(min)"] - dataFameC[j, "Average Rt(min)"])
      
      i <- i + 1
    }
    
    # 数据整理
    a <- ncol(data3)
    data3 <- data3[, c(1:2, a, 3:(a - 1))]
    savexlsx1(filename = savefile, sheet = "RI数据矩阵", data = data3)
  } else {
    print("未进行RI计算")
  }
  
  colnames(data3)[colnames(data3) == "Alignment ID"] <- "ID"
  colnames(data3)[colnames(data3) == "Average Rt(min)"] <- "Retention time (min)"
  colnames(data3)[colnames(data3) == "Quant mass"] <- "m/z"
  colnames(data3)[colnames(data3) == "Total score"] <- "Score"
  
  # GCMS level 
  data3[,"level"] <- "Level 3"
  data3[(data3[,"Score"] > 80),"level"] <- "Level 2"
  data3[(data3[,"Score"] > 85),"level"] <- "Level 1"
  
  # info整理
  nameorder <- c("ID","Retention time (min)","Average RI","m/z",
                 "Metabolites","Score","cid","CAS","HMDB","METLIN",
                 "Lipidmaps","KEGG","PubChem","ChEBI","smiles","InChIKey",
                 "Formula","Super Class","Class","Sub Class","level")
  data3 <- data3[unique(c(intersect(nameorder,names(data3)),names(data3)[!names(data3) %in% nameorder]))]
  
  # 按order排序
  data3 <- data3[order(data3$Score, decreasing = T),]
  data3 <- data3[order(data3$level, decreasing = F),]

  savexlsx1(filename = savefile, sheet = "数据矩阵", data = data3)
  return(data3)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",default = "Area.txt", help = "GCMS搜库下机数据")
  parser$add_argument("-m","--mode",default = "fame", help = "归一化方式,area(峰面积)、fame(内标分段归一化)、single(单一内标归一化)")
  parser$add_argument("-p","--projectclass",default = "gc", help = "数据类型,spme(顶空)、gc(非靶)、wax(蜡质)")
  parser$add_argument("-fr","--famersd",default = 0.1, type = "double",help = "内标分段归一化的内标rsd筛选标准")
  parser$add_argument("-fa","--fameaverage",default = F,help = "是否对内标中0值进行均值填充",action='store_true')
  parser$add_argument("-s","--score",default = 70, type = "double",help = "定性打分筛选标准")
  parser$add_argument("-pl","--peaklimit",default = 5000, type = "double",help = "峰强最小限制")
  parser$add_argument("-sa","--savefile",default = "数据矩阵.xlsx",help = "文件保存名称")
  parser$add_argument("-r","--RI",default = F,help = "是否进行RI运算",action='store_true')
  parser$add_argument("-sf","--singlefame",default = NULL,help = "单内标名称")
  parser$add_argument("-dl","--deallist",default = NULL,help = "需处理的物质名称",nargs = "+")
  
  args <- parser$parse_args()
  
  result <- do.call(what = datafetch_untargetgcms,args = args) 
  
}
