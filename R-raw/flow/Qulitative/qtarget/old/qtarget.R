#!/opt/conda/bin/Rscript

# utils function
#' @export
colnames_refine <- function (colnames) {
  res <- vector()
  for (colname in colnames){
    arr <- strsplit(colname,"-",fixed = T)[[1]]
    if (length(arr) <= 2){
      res <- c(res, colname)
    } else {
      suffix <- arr[3:length(arr)]
      new_colname <- paste(suffix,collapse = "-")
      res <- c(res, new_colname)
    }
  }
  return(res)
}

# 1.数据读取----------------------------------------------------------
#' @export
getmzmlData2 <- function(arg,
                         methodfrom = "^animal",
                         methodpath =  databasepath,
                         ...) {
  
  data <- list()
  mzmlData <- list()
  method <- list()
  
  arg1 <- arg
  
  for (mode1 in unique(unlist(arg1["模式"]))) {
    if (mode1 == "neg") {
      methodFile <- list.files(
        path = methodpath,
        pattern = paste0(methodfrom, "-neg"), full.names = T
      )
      methodFile <- methodFile[length(methodFile)]
    } else if (mode1 == "pos") {
      methodFile <- list.files(
        path = methodpath,
        pattern = paste0(methodfrom, "-pos"), full.names = T
      )
      methodFile <- methodFile[length(methodFile)]
    } else {
      stop('path must be "neg" or "pos"')
    }
    
    database <- read.csv(methodFile, header = T, sep = "\t", check.names = F)
    database1 <- database %>%
      mutate(pre_pro = paste0(`Precursor (Q1) Mass (Da)`, "_", `Fragment (Q3) Mass (Da)`))
    print(paste0("mzml文件读取中:", mode1))
    mzml <- unlist(arg1[which(arg1["模式"] == mode1), c("mzml文件名称")])
    pd <- arg1[which(arg1["模式"] == mode1), c("组别", "mzml文件名称")]
    raw_data <- MSnbase::readSRMData(files = mzml, pdata = pd)
    # 用伪数据 fill na rt
    for (nrow in 1:dim(raw_data)[1]){
      for (ncol in 1:dim(raw_data)[2]){
        if (length(raw_data[nrow,ncol]) <= 1){
          raw_data[nrow,ncol]@rtime <- c(0.11,0.12,0.13)
          raw_data[nrow,ncol]@intensity <- c(0,0,0)
        }
      }
    }
    mzmlData <- c(mzmlData, list(raw_data))
    method <- c(method, list(database1))
  }
  
  names(mzmlData) <- unique(unlist(arg1["模式"]))
  names(method) <- unique(unlist(arg1["模式"]))
  data <- c(list(arg1), list(mzmlData), list(method))
  names(data) <- c("args", "mzmlData", "method")
  return(data)
}

# 2.峰检测-------------------------------------------
# centwave算法,关键参数：
# 1.peakwidth 色谱峰宽的预期范围
# 2.ppm 对应于一个色谱峰的质心 m/z 值的最大预期偏差；这通常远大于制造商指定的 ppm
# 注意：SRM中只有expectRT和特定离子对，即，只能根据rt来计算峰面积
#' @export
getPeaks2 <- function(mzmlData,
                      peakWidth = c(0.2, 0.4),
                      noise = 2000,
                      snthresh = 3) {
  res <- list()
  snthresh <- ceiling(snthresh)
  print(paste0("噪声为：", noise))
  print(paste0("信噪比为：", snthresh))
  cwp <- CentWaveParam(
    peakwidth = peakWidth,
    noise = noise,
    snthresh = snthresh,
    fitgauss = T
  )
  
  for (name1 in names(mzmlData[["mzmlData"]])) {
    raw_data <- mzmlData[["mzmlData"]][[name1]]
    print(paste0("色谱峰提取中：", name1))
    
    xdata <- findChromPeaks(raw_data, cwp, BPPARAM = bpparam())
    xdata1 <- chromPeaks(xdata)
    res1 <- data.frame(xdata1)[c("row", "column", "rt", "sn", "into")]
    res1[is.na(res1)] <- 0
    res1["into"] <- res1["into"] * 60
    res1 <- res1 %>%
      mutate(row_col = paste0(row, "_", column))
    res <- c(res, list(res1))
  }
  names(res) <- names(mzmlData[["mzmlData"]])
  mzmlData[["peaks"]] <- res
  return(mzmlData)
}

# 3.数据处理------------------------------------------------------
# 根据row、column信息获取pre\pro\name\expectRT信息
# q1_q3和rt做数据筛选
#' @export
peakProcess2 <- function(mzmlData) {
  peakInfo <- list()
  for (name1 in names(mzmlData[["peaks"]])) {
    
    res <- mzmlData[["peaks"]][[name1]]
    raw_data <- mzmlData[["mzmlData"]][[name1]]
    database1 <- mzmlData[["method"]][[name1]]
    
    for (i in 1:dim(res)[1]) {
      row1 <- res[i, "row"]
      col1 <- res[i, "column"]
      pre <- raw_data[row1, col1]@precursorMz[1]
      pre1 <- round(pre, 2)
      pro <- raw_data[row1, col1]@productMz[1]
      pro1 <- round(pro, 2)
      
      pre_pro <- paste0(pre1, "_", pro1)
      rt <- res[i, "rt"]
      hit_record <- database1[database1["pre_pro"] == pre_pro,]
      
      if (dim(hit_record)[1] == 1) {
        expectRT <- hit_record[1, "Retention Time (min)"]
        name <- hit_record[1, "Name"]
        
        count1 <- which(abs(rt - expectRT) < 0.3)
        if (length(count1) == 0) {
          res[i, "pre_pro"] <- pre_pro
          res[i, "rtDrift"] <- "drift"
        } else if (length(count1) == 1) {
          res[i, "pre_pro"] <- pre_pro
          res[i, "name"] <- name
          res[i, "expectRT"] <- expectRT
          res[i, "rtDrift"] <- "noDrift" }
      } else if (dim(hit_record)[1] > 1) {
        warning("Records with same preMZ and proMZ in database were found!")
        expectRT <- hit_record[, "Retention Time (min)"]
        name <- hit_record[, "Name"]
        
        count1 <- which(abs(rt - expectRT) < 0.3)
        if (length(count1) == 0) {
          res[i, "pre_pro"] <- pre_pro
          res[i, "rtDrift"] <- "drift"
        } else {
          res[i, "pre_pro"] <- pre_pro
          res[i, "name"] <- name[count1[1]]
          res[i, "expectRT"] <- expectRT[count1[1]]
          res[i, "rtDrift"] <- "sameNameNearbyExpectRT"}
      } else {
        res[i, "pre_pro"] <- pre_pro
        res[i, "rtDrift"] <- "drift"
      }
    }
    peakInfo <- c(peakInfo, list(res))
  }
  names(peakInfo) <- names(mzmlData[["peaks"]])
  mzmlData[["peakInfo"]] <- peakInfo
  return(mzmlData)
}


# 4.数据合并-----------------------------------
#' @export
dataCombine2 <- function(mzmlData1,
                         mzmlData2) {
  combined <- list()
  for (name1 in names(mzmlData1[["peakInfo"]])) {
    temp1 <- mzmlData1[["peakInfo"]][[name1]]
    temp2 <- temp1[!is.na(temp1["name"]), ]
    temp2["level"] <- "I"
    mzml_cols <- max(unique(temp2["column"]))
    for (row1 in unlist(unique(temp2["row"]))) {
      temp <- temp2[temp2["row"] == row1, ]
      ratio <- dim(temp)[1] / mzml_cols
      temp2[temp2["row"] == row1, "naRatio"] <- 1 - ratio
    }
    temp3 <- mzmlData2[["peakInfo"]][[name1]]
    temp4 <- temp3[!is.na(temp3["name"]), ]
    temp4["level"] <- "II"
    mzml_cols <- max(unique(temp2["column"]))
    for (row1 in unlist(unique(temp4["row"]))) {
      temp <- temp4[temp4["row"] == row1, ]
      ratio <- dim(temp)[1] / mzml_cols
      temp4[temp4["row"] == row1, "naRatio"] <- 1 - ratio
    }
    temp <- rbind(temp2, temp4)
    
    res <- temp %>%
      mutate(row_col = paste0(row, "_", column)) %>%
      arrange(naRatio, level) %>%
      distinct(row_col, .keep_all = T) %>%
      mutate(name_col = paste0(name, "_", column)) %>%
      arrange(desc(sn)) %>%
      distinct(name_col, .keep_all = T)
    res1 <- res[res["sn"] >= 2, ]
    res1 <- res1[res1["naRatio"] < 0.8, ]
    combined <- c(combined, list(res1))
  }
  names(combined) <- names(mzmlData1[["peakInfo"]])
  mzmlData1[["combined"]] <- combined
  return(mzmlData1)
}

#' @export
getRes2 <- function(mzmlData,outputpath) {
  peakRes <- list()
  temp <- mzmlData[["combined"]]
  arg1 <- mzmlData[["args"]]
  for (name1 in names(temp)) {
    d1 <- temp[[name1]]
    colnames(d1)
    temp1 <- d1 %>%
      arrange(column) %>%
      pivot_wider(id_cols = "name", names_from = "column", values_from = "into")
    temp_cols <- c("name", 1:(dim(temp1)[2] - 1))
    temp1 <- temp1[temp_cols]
    temp2 <- mzmlData[["method"]][[name1]][c(1, 2, 3, 4)]
    names(temp2) <- c("Identified metabolite", "Precursor m/z", "Product m/z", "Retention time")
    res <- merge(temp2, temp1, by.x = "Identified metabolite", by.y = "name", all.y = T)
    colnames(res) <- c(colnames(temp2), unlist(arg1[arg1[1] == name1, 3]))
    res[is.na(res)] <- 0
    res["sum"] <- rowSums(res[-c(1:4)])
    # colnames(res)
    res1 <- res %>%
      arrange(desc(sum))
    # len <- dim(res1)[1] * 0.95
    len <- dim(res1)[1]
    sum_col <- dim(res1)[2]
    res1 <- res1[1:len, -sum_col]
    colnames(res1)<-colnames_refine(colnames(res1))
    peakRes <- c(peakRes, list(res1))
    res2 <- as.data.frame(t(res1))
    write.table(res2, file = paste0(outputpath, "/", name1, ".txt"), quote = F, sep = "\t", dec = ".", row.names = T, col.names = F)
  }
  names(peakRes) <- names(temp)
  mzmlData[["peakRes"]] <- peakRes
  openxlsx::write.xlsx(mzmlData[["peakRes"]], paste0(outputpath,"/", "res.xlsx"), overwrite = T)
  return(mzmlData)
}

# 6.绘图
#' @export
drawPic <- function(mzmlData,outputpath) {
  library(RColorBrewer)
  for (mode1 in names(mzmlData[["combined"]])) {
    print(paste0("开始绘图：",mode1))
    temp <- mzmlData[["combined"]][[mode1]]
    dir_fig <- paste0(outputpath,"/",mode1, "_有峰信息")
    unlink(dir_fig, force = T, recursive = T)
    dir.create(dir_fig)
    
    dir_fig1 <- paste0(outputpath,"/",mode1, "_无峰信息")
    unlink(dir_fig1, force = T, recursive = T)
    dir.create(dir_fig1)
    mz_groups <- unlist(unique(mzmlData[["args"]][mzmlData[["args"]]["模式"]==mode1,2]))
    mz_color_set <- brewer.pal(length(mz_groups), "Set1")
    mz_colors <- c()
    for (g in 1:length(mz_groups)) {
      temp1 <- rep(mz_color_set[g], sum(mzmlData[["args"]][mzmlData[["args"]]["模式"]==mode1,2] == mz_groups[g]))
      mz_colors <- c(mz_colors, temp1)
    }
    for (row1 in 1:dim(mzmlData[["mzmlData"]][[mode1]])[1]) {
      if (any(temp["row"] == row1)) {
        name <- unique(temp[temp["row"] == row1, "name"])
        pre_pro <- unique(temp[temp["row"] == row1, "pre_pro"])
        expectRT <- unique(temp[temp["row"] == row1, "expectRT"])
        rt1 <- round(mean(temp[temp["row"] == row1, "rt"]), 2)
        title1 <- paste0(name, "_", pre_pro, "_", expectRT, "_", rt1)
        jpeg(paste0(dir_fig, "/", row1, ".jpg"))
        plot(mzmlData[["mzmlData"]][[mode1]][row1], col = mz_colors, main = "")
        legend("topleft", mz_groups, col = mz_color_set, lty = c(1, 1, 1))
        title(title1)
        dev.off()
      } else {
        pre <- mzmlData[["mzmlData"]][[mode1]][row1, 1]@precursorMz[1]
        pro <- mzmlData[["mzmlData"]][[mode1]][row1, 1]@productMz[1]
        pre_pro <- paste0(pre, "_", pro)
        # 多个匹配报错
        databaseRT <- mzmlData[["method"]][[mode1]][mzmlData[["method"]][[mode1]]["pre_pro"] == pre_pro, "Retention Time (min)"][1]
        jpeg(paste0(dir_fig1, "/", row1, ".jpg"))
        plot(mzmlData[["mzmlData"]][[mode1]][row1], col = mz_colors, main = "")
        legend("topleft", mz_groups, col = mz_color_set, lty = c(1, 1, 1))
        title(paste0(pre_pro, "_", databaseRT))
        dev.off()
      }
    }
  }
}

#' @export
QtargetQualitative2 <- function(mzmlpath = "./mzml",
                                sampleinfo = 'sampleinfo.xlsx',
                                methodfrom = "^animal",
                                methodpath = "./db",
                                outputpath = "./output",
                                ...) {
  library(xcms)
  library(mzR)
  library(tidyverse)
  library(magrittr)
  library(readxl)
  # if (file.exists(outputpath)){
  #   warning(paste0("文件夹",outputpath,"已存在!"))
  # } else {
  #   dir.create(outputpath)
  # }
  data <- readxl::read_xlsx(sampleinfo)
  arg <- data %>% mutate(`mzml文件名称` = paste0(mzmlpath,"/",`模式`,"/",`mzml文件名称`))
  
  arg <- arg[!is.na(arg$`组名`), ]
  arg <- arg[!is.na(arg$`mzml文件名称`), ]
  
  if (nrow(arg) == 0) {
    stop("无定性数据")
  }
  
  mzmlData <- getmzmlData2(arg = arg,
                           methodfrom = methodfrom,
                           methodpath = methodpath)
  
  mzmlData1 <- getPeaks2(mzmlData = mzmlData,
                         peakWidth = c(0.2, 0.4), # 峰宽
                         noise = 2000, # 噪声
                         snthresh = 3) # 信噪比
  mzmlData1 <- peakProcess2(mzmlData = mzmlData1)
  
  # 放宽条件
  mzmlData2 <- getPeaks2(mzmlData = mzmlData,
                         peakWidth = c(0.2, 0.4),
                         noise = 500,
                         snthresh = 2)
  
  mzmlData2 <- peakProcess2(mzmlData = mzmlData2)
  
  mzmlData <- dataCombine2(mzmlData1, mzmlData2)
  
  mzmlData <- getRes2(mzmlData,outputpath)
  
  drawPic(mzmlData,outputpath)
  
  negfile <- list.files(path = paste0(outputpath,"/"), pattern = "^neg.txt$")
  if (length(negfile) == 0) {
    warning("负离子模式拟靶向结果未查询到")
  }
  posfile <- list.files(path = paste0(outputpath,"/"), pattern = "^pos.txt$")
  if (length(posfile) == 0) {
    warning("正离子模式拟靶向结果未查询到")
  }
  if (length(negfile) == 0 & length(posfile) == 0) {
    stop("正负离子模式拟靶向结果均未查询到，暂停分析，请查找问题原因")
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-m","--mzmlpath",default = "/data/test_data/target/test-project", help = "拟靶向项目存放mzml文件的pos和neg文件夹的父目录路径,输入路径末尾不要加 '/'")
  parser$add_argument("-s","--sampleinfo",default = "sampleinfo.xlsx", help = "拟靶向项目样本基本信息xlsx,sampleinfo.xlsx必须包含模式、组别、组名、mzml文件名称字段")
  parser$add_argument("-f","--methodfrom",default = "^animal", help = "数据库文件开头的字符，^animal")
  parser$add_argument("-p","--methodpath",default = "/data/hstore4/database/qualitative/QtargetLipid", help = "拟靶向项目数据库的存放目录,default:/data/hstore4/database/qualitative/QtargetLipid输入路径末尾不要加 '/'")
  parser$add_argument("-o","--outputpath",default = "/output", help = "拟靶向项目结果的输出目录,default[./output]输入路径末尾不要加 '/'")
  
  args <- parser$parse_args()
  
  result <- do.call(what = QtargetQualitative2,args = args)
}
