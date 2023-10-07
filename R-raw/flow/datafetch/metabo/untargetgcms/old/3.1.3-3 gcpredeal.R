#' @export
gcpredeal <- function(filename = NULL,
                      mode = "fame", # area峰面积,fame内标分段归一化,single单一内标归一化
                      projectclass = "gc", # spme(顶空) or gc(非靶gc) or wax(蜡质)
                      famersd = 0.1, # 内标rsd筛选标准0.15
                      score = 70,
                      peaklimit = 5000,
                      savefile = "数据矩阵.xlsx",
                      savesheet = "数据矩阵",
                      fameaverage = F,
                      singlefame = "Decanoic Acid, Ethyl Ester", # 单内标名字Decanoic acid, ethyl ester葵酸乙酯
                      deallist = c("Fame C", "Unknown", "Tbs Compound", "Z [^a-z]", "3,4-Dichlorophenylalanine",
                                   "fluo", "chlo", "brom", "iod", "sulf", "silane","sil"),
                      spmedeallist = c("tms derivative", "Pyridine", "Propoxur", "Ethanol", "tbdms derivative"),
                      waxdeallist = c("pyridine", "propoxur", "ethanol", "tbdms derivative"),
                      classfile = "/data/hstore4/database/gcpredeal/LUG-V3.3-HMDB-CLASS.xlsx",
                      RI = T ,# 是否进行RI计算
                      ...) {
  
  datainfo1 <- data.frame(filename, mode, projectclass, famersd, score,peaklimit)

  if (projectclass == "spme") {
    datainfo2 <- data.frame(c(deallist, spmedeallist))
  }else if ( projectclass == "gc"){
    datainfo2 <- data.frame(deallist)
  }else if ( projectclass == "wax"){
    datainfo2 <- data.frame(c(deallist, waxdeallist))
  }

  colnames(datainfo2) <- "deallist"
  savexlsx1(filename = "GC预处理info.xlsx", sheet = "info1", data = datainfo1)
  savexlsx1(filename = "GC预处理info.xlsx", sheet = "info2", data = datainfo2)
  rm(datainfo1, datainfo2)

  # 数据读入
  spmedata <- readdata(filename = filename, rowNames = F, header = F, skip = 0)
  dataclass <- readdata(filename = classfile, sheet = 1)
  dataclass$Name <- stringr::str_to_title(dataclass$Name)
  colnames(dataclass) <- c("Metabolites", "HMDB", "Super Class", "Class", "Sub Class")
  famenorm <- data.frame()

  # 数据整理
  colnames(spmedata) <- spmedata[5, ]
  colnames(spmedata)[5] <- "Metabolites"
  #首字母大写
  spmedata$Metabolites <- stringr::str_to_title(spmedata$Metabolites)
  singlefame <- stringr::str_to_title(singlefame)
  # 删除空的表头
  dataraw <- spmedata[6:(length(spmedata[, 1])), 1:(max(max(which(spmedata[2, ] == "Sample"), max(which(spmedata[2, ] == "QC")))))]
  rownames(dataraw) <- dataraw[, 1]
  # 定性数据提出
  dataID <- dataraw[, c(1, 2, 4, 5, 11, 19)]
  dataID[, c(1, 2, 3, 6)] <- as.numeric(unlist(dataID[, c(1:3, 6)]))
  # 定量数据提出
  dataM <- data.frame(lapply(dataraw[, 29:length(dataraw[1, ])], as.numeric), check.names = F)
  #低于5000的峰归0
  dataM[dataM<= peaklimit] <- 0
  # 原始数据保存
  dataraw <- data.frame(dataID, dataM, check.names = F)
  savexlsx1(filename = savefile, sheet = "原始数据矩阵", data = dataraw)
  # 内标数据提取
  dataFameC <- dataraw[grep("Fame C", dataraw$Metabolites), ]
  singlefame_area <- dataraw[grep(singlefame, dataraw$Metabolites), ]
  singlefame_area <- arrange(singlefame_area, desc(`Total score`))
  dataFameCM <- dataFameC[, 7:length(dataFameC)]
  singlefameM <- singlefame_area[!duplicated(singlefame_area$Metabolites), 7:length(singlefame_area)]

  # 内标数据0值检查
  if (mode == "fame") {
    dataFameCM[dataFameCM == 0] <- NA
    n <- sum(is.na(dataFameCM))
    if (n != 0) {
      print(paste0("Fame内标异常:", n, "个内标为空值,请检查内标!!!检查无误后请使用fameaverage = T,后续将以该内标所有样本均值替代空值"))
    } else {
      print("内标无异常,请使用fameaverage = F继续运行")
    }
  }else if (mode == "single") {
    singlefameM[singlefameM == 0] <- NA
    n <- sum(is.na(singlefameM))
    if (n != 0) {
      print(paste0("single内标异常:", n, "个内标为空值,请检查内标!!!检查无误后请使用fameaverage = T,后续将以该内标所有样本均值替代空值"))
    } else {
      print("内标无异常,请使用fameaverage = F继续运行")
    }
  }

  # 内标数据0值处理
  if (fameaverage) {
    if (mode == "fame") {
      for (i in 1:length(dataFameCM[, 1])) {
        dataFameCM[i, c(which(is.na(dataFameCM[i, ])))] <- apply(dataFameCM[i, -c(which(is.na(dataFameCM[i, ])))], 1, mean)
      }
    } else if (mode == "single") {
      for (i in 1:length(singlefameM[, 1])) {
        singlefameM[i, c(which(is.na(singlefameM[i, ])))] <- apply(singlefameM[i, -c(which(is.na(singlefameM[i, ])))], 1, mean)
      }
    }
  }else {
    print("fameaverage = F,不进行内标处理")
  }

  # 3种归一化
  if (mode == "area") {
    print("当前正在使用的归一化方法:峰面积归一化")
    # 计算总峰面积
    allarea <- apply(X = dataM, MARGIN = 2, FUN = sum)
    # 归一化后ID和M文件合并并保存
    data1 <- data.frame(dataID, t(t(dataM) / allarea * 10000), check.names = F)
    savexlsx1(filename = savefile, sheet = "峰面积归一化", data = data1)
  }else if (mode == "fame") {
    if (sum(is.na(dataFameCM)) > 0) {
      return(message("内标有问题,进程终止,请检查内标"))
    }

    else {
      # 内标分段归一化
      print("当前正在使用的归一化方法:内标分段归一化")
      # 计算所有样本内标rsd并保存
      datarsd <- data.frame(
        `mean-ALL` = apply(dataFameC[, 7:length(dataFameC)], MARGIN = 1, FUN = mean),
        `sd-ALL` = apply(dataFameC[, 7:length(dataFameC)], MARGIN = 1, FUN = sd),
        check.names = F
      )
      # 计算qc样本的内标rsd并保存
      # datarsd <- data.frame(`mean-qc` = apply(dataFameC[,grep("QC",colnames(dataFameC),ignore.case = T)], MARGIN = 1, FUN = mean),
      #                      `sd-qc` = apply(dataFameC[,grep("QC",colnames(dataFameC),ignore.case = T)], MARGIN = 1, FUN = sd),
      #                      check.names = F)
      # dataFameC$'rsd-qc' <- datarsd$`sd-qc`/datarsd$`mean-qc`
      dataFameC$"rsd-ALL" <- datarsd$`sd-ALL` / datarsd$`mean-ALL`
      dataFameC <- dataFameC[, c(1:6, length(dataFameC), 7:(length(dataFameC) - 1))]
      savexlsx1(filename = savefile, sheet = "FAME", data = dataFameC)

      # 内标分段归一化
      # 内标rsd筛选
      datafamersd <- dataFameC %>% filter(dataFameC$`rsd-ALL` < famersd)
      datafameM <- datafamersd[8:length(datafamersd[1, ])]
      timelist <- c(0, datafamersd$`Average Rt(min)`)
      # 分段归一(part1)
      for (i in 2:(length(timelist))) {
        datasub <- dataM %>% filter(dataraw$`Average Rt(min)` > timelist[i - 1], dataraw$`Average Rt(min)` <= timelist[i])
        famenorm <- rbind(famenorm, t(t(datasub) / as.numeric(datafameM[i - 1, ])))
      }
      # 分段归一(part2)
      datasub <- dataM %>% filter(dataraw$`Average Rt(min)` > timelist[length(timelist)])
      famenorm <- rbind(famenorm, datasub / as.numeric(datafameM[length(datafameM[, 1]), ]))
      # 内标分段归一化结果M与ID合并
      data1 <- data.frame(dataID, famenorm, check.names = F)
      savexlsx1(filename = savefile, sheet = "内标分段归一化", data = data1)


    }
  } else if (mode == "single") {
    if (sum(is.na(singlefameM)) > 0) {
      return(message("内标有问题,进程终止,请检查内标"))
    }
    else {
      # 单一内标归一化
      print("当前正在使用的归一化方法:单一内标归一化")
      data1 <- data.frame(dataID, t(t(dataM) / as.numeric(singlefameM)), check.names = F)
      savexlsx1(filename = savefile, sheet = paste0(singlefame, "归一化"), data = data1)
    }
  }


  # 筛选
  # 1.编号处理
  data1$Metabolites <- gsub(";[0-9]+$", "", data1$Metabolites)
  # 2.删除unknown/z /tbs等代谢物
  data1$Metabolites <- tolower(data1$Metabolites)
  deallist <- tolower(deallist)
  singlefame <- tolower(singlefame)
  spmedeallist <- tolower(spmedeallist)
  waxdeallist <- tolower(waxdeallist)
  if (projectclass == "spme") {
    data2 <- data1[-c(as.numeric(grep(paste0(c(deallist, spmedeallist, singlefame), collapse = "|"), data1$Metabolites))), ]
  }else if ( projectclass == "gc"){
    data2 <- data1[-c(as.numeric(grep(paste0(c(deallist, singlefame), collapse = "|"), data1$Metabolites))), ]
  } else if ( projectclass == "wax"){
    data2 <- data1[-c(as.numeric(grep(paste0(c(waxdeallist,deallist, singlefame), collapse = "|"), data1$Metabolites))), ]
  }
  data2$Metabolites <- stringr::str_to_title(data2$Metabolites)
  # 3.删除打分低于score代谢物
  data2 <- data2 %>% filter(data2$`Total score` > score)
  # 4.1去除名字重复
  data2 <- arrange(data2, desc(`Total score`))
  data2 <- data2[!duplicated(data2$Metabolites), ]
  # 4.2去除inchikey重复
  data2 <- data2[!duplicated(data2$INCHIKEY), ]
  # 5去除inchikey列
  data2 <- data2[, -5]
  if (projectclass %in% c("spme","wax")) {
    savexlsx1(filename = savefile, sheet = savesheet, data = data2)
    data3 <- data2
  } else {
    data3 <- merge(data2, dataclass, by = "Metabolites", all.x = T)
    data3 <- data3[, c(2:4, 1, 5, (length(data3) - 3):length(data3), 6:(length(data3) - 4))]
    savexlsx1(filename = savefile, sheet = savesheet, data = data3)
  }
  finaldata <- data3
  # RI指数计算
  if (RI) {
    datari <- readdata(filename = "/data/hstore4/database/gcpredeal/RI.xlsx", sheet = 1)
    dataFameC <- merge(x = dataFameC, y = datari, by = "Metabolites")
    data3 <- data3[order(data3$`Average Rt(min)`), ]
    i <- 1
    while (i <= dim(data3)[1]) {
      data3[i, "Average Rt(min)"]

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
    finaldata <- data3[, c(1:2, a, 3:(a - 1))]
    savexlsx1(filename = savefile, sheet = "RI数据矩阵", data = finaldata)
  } else {
    print("未进行RI计算")
  }
  return(finaldata)
}
