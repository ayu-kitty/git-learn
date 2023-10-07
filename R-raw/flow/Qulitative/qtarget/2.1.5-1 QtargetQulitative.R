#!/opt/conda/bin/Rscript

# 测试
# library("lmbio")
# setwd("/data/hstore1/database/test/2023-03-14拟靶向定性")
# filename <- list.files(path = "neg",full.names = T)
# mode <- "neg"
# methodpath <- selectfile(path = "database/qualitative/QtargetLipid",file = paste0("qualitative-",mode,".txt"))
# outputpath <- "../raw/搜库数据"

#' 拟靶数据定性
#'
#' @param filename 数据名
#' @param methodpath 方法路径
#' @param outputpath 输出路径
#' @param ...
#'
#' @export
QtargetQualitative <- function(filename,
                               methodpath = selectfile(path = "database/qualitative/QtargetLipid",
                                                       file = paste0("qualitative-",mode,".txt")),
                               outputpath = "../raw/搜库数据",
                               mode = "neg",
                               ...){
  print("拟靶向定性流程-2023-03-14")
  
  mzmlData <- getmzmlData(filename = filename,
                          methodpath = methodpath) 
  
  mzmlData <- getmzmlPeaks(mzmlData = mzmlData,
                           peakWidth = c(0.08, 0.4),
                           noise = 500,
                           snthresh = 2)
  
  mzmlData <- peakProcess(mzmlData = mzmlData)
  
  mzmlData <- getRes(mzmlData = mzmlData)
  
  filepath <- paste0(outputpath,"/",mode,".txt")
  
  savetxt(data = mzmlData$peakRes,filename = filepath,
          row.names = T,col.names = F)
  
  mzmlData[["ouputpath"]] <- filepath
  
  return(mzmlData)
}


# 1.数据读取----------------------------------------------------------
getmzmlData <- function(filename,
                        methodpath) {
  suppressMessages(library("xcms"))
  suppressMessages(library("mzR"))
  suppressMessages(library("tidyverse"))
  suppressMessages(library("magrittr"))
  
  database <- readdata(methodpath)
  database1 <- database %>%
    mutate(pre_pro = paste0(`Precursor (Q1) Mass (Da)`, "_", `Fragment (Q3) Mass (Da)`))
  
  raw_data <- readdata(filename = filename)
  
  # 用伪数据 fill na rt
  for (nrow in 1:dim(raw_data)[1]){
    for (ncol in 1:dim(raw_data)[2]){
      if (length(raw_data[nrow,ncol]) <= 1){
        raw_data[nrow,ncol]@rtime <- c(0.11,0.12,0.13)
        raw_data[nrow,ncol]@intensity <- c(0,0,0)
      }
    }
  }
  
  data <- list(filename = filename,
               methodpath = methodpath,
               method = database1,
               mzmlData = raw_data)
  
  return(data)
}

# 2.峰检测-------------------------------------------
# centwave算法,关键参数：
# 1.peakwidth 色谱峰宽的预期范围
# 2.ppm 对应于一个色谱峰的质心 m/z 值的最大预期偏差；这通常远大于制造商指定的 ppm
# 注意：SRM中只有expectRT和特定离子对，即，只能根据rt来计算峰面积
getmzmlPeaks <- function(mzmlData,
                         peakWidth = c(0.2, 0.4),
                         noise = 2000,
                         snthresh = 3) {
  
  snthresh <- ceiling(snthresh)
  print(paste0("噪声为:", noise))
  print(paste0("信噪比为:", snthresh))
  cwp <- CentWaveParam(peakwidth = peakWidth,
                       noise = noise,
                       snthresh = snthresh,
                       firstBaselineCheck = F,
                       fitgauss = F)
  
  raw_data <- mzmlData[["mzmlData"]]
  
  xdata <- findChromPeaks(raw_data, cwp, BPPARAM = bpparam())
  xdata1 <- chromPeaks(xdata)
  res1 <- data.frame(xdata1)[c("row", "column", "rt", "sn", "into")]
  res1[is.na(res1)] <- 0
  res1["into"] <- res1["into"] * 60
  res1 <- res1 %>%
    mutate(row_col = paste0(row, "_", column))
  
  res1 <- res1[res1$sn >= snthresh,]
  res1 <- res1[res1$into >= noise,]
  
  mzmlData[["peaks"]] <- res1
  
  return(mzmlData)
}

# 3.数据处理------------------------------------------------------
# 根据row、column信息获取pre\pro\name\expectRT信息
# q1_q3和rt做数据筛选
peakProcess <- function(mzmlData) {
  
  res <- mzmlData[["peaks"]]
  raw_data <- mzmlData[["mzmlData"]]
  database1 <- mzmlData[["method"]]
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
    }else if(dim(hit_record)[1] == 0){
      warning("数据库未找到")
      next
    }else{
      # warning("数据库有重复")
      expectRT <- hit_record[1, "Retention Time (min)"]
      name <- hit_record[1, "Name"]
    }
    
    count1 <- which(abs(rt - expectRT) < 0.5)
    if (length(count1) == 0) {
      res[i, "pre_pro"] <- pre_pro
      res[i, "rtDrift"] <- "Drift"
    } else if (length(count1) == 1) {
      res[i, "pre_pro"] <- pre_pro
      res[i, "name"] <- name
      res[i, "expectRT"] <- expectRT
      res[i, "rtDrift"] <- "noDrift" 
    } 
  }
  
  res1 <- res[res$rtDrift == "noDrift",]
  res2 <-  res1[grepl(pattern = "IS$",x = res1$name),]
  if(dim(res2)[1] > 5){
    res1 <- res2
  }
  
  for ( i in 1:dim(res1)[1]) {
    res2 <- res1[res1$name == res1[i,"name"],]
    res2 <- res2[abs(res2$rt - res1[i,"rt"]) < 0.1,]
    res2 <- res2[!duplicated(res2$column),]
    res1[i,"num"] <- dim(res2)[1]
  }
  res2 <- res1
  res2 <- res2[order(res2$into,decreasing = T),]
  res2 <- res2[order(res2$num,decreasing = T),]
  res2 <- res2[!duplicated(res2$name),]
  x <- res2$expectRT
  y <- res2$rt
  lmdata <- lm(y~x,weights = 1/x)
  database1[,"newtime"] <- predict(lmdata,data.frame(x=c(database1$`Retention Time (min)`)))
  database1[,"Retention Time (min)"] <- predict(lmdata,data.frame(x=c(database1$`Retention Time (min)`)))
  
  res <- mzmlData[["peaks"]]
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
    }else if(dim(hit_record)[1] == 0){
      warning("数据库未找到")
      next
    }else{
      # warning("数据库有重复")
      expectRT <- hit_record[1, "Retention Time (min)"]
      name <- hit_record[1, "Name"]
    }
    
    count1 <- which(abs(rt - expectRT) < 0.3)
    if (length(count1) == 0) {
      res[i, "pre_pro"] <- pre_pro
      res[i, "rtDrift"] <- "Drift"
    } else if (length(count1) == 1) {
      res[i, "pre_pro"] <- pre_pro
      res[i, "name"] <- name
      res[i, "expectRT"] <- expectRT
      res[i, "rtDrift"] <- "noDrift" 
    } 
  }
  
  res1 <- res[res$rtDrift == "noDrift",]
  
  for ( i in 1:dim(res1)[1]) {
    res2 <- res1[res1$name == res1[i,"name"],]
    res2 <- res2[abs(res2$rt - res1[i,"rt"]) < 0.1,]
    res2 <- res2[!duplicated(res2$column),]
    res1[i,"num"] <- dim(res2)[1]
  }
  res1[,"diffrt"] <- abs(res1$rt-res1$expectRT)
  
  res2 <- res1
  res2 <- res2[order(res2$into,decreasing = T),]
  res2 <- res2[order(res2$diffrt),]
  res2 <- res2[order(res2$num,decreasing = T),]
  res2 <- res2[!duplicated(res2$name),]
  
  res3 <- NULL
  
  for ( i in 1:dim(res2)[1]) {
    res4 <- res[res$row == res2[i,"row"],]
    res4 <- res4[abs(res4$rt - res2[i,"rt"]) < 0.15,]
    res4 <- res4[order(abs(res4$rt - res2[i,"rt"])),]
    res4 <- res4[!duplicated(res4$column),]
    res4$name <- res2[i,"name"]
    res4$expectRT <- res2[i,"rt"]
    res4$rtDrift <- "noDrift"
    res3 <- rbind(res3,res4)
  }
  res3 <- res3[order(res3$column),]
  res3 <- res3[order(res3$row),]
  
  mzmlData[["peakInfo"]] <- res3
  return(mzmlData)
}

getRes <- function(mzmlData) {
  
  d1 <- mzmlData[["peakInfo"]]
  temp1 <- d1 %>%
    arrange(column) %>%
    pivot_wider(id_cols = "name", names_from = "column", values_from = "into")
  colnames(temp1) <- c("name", colnames(mzmlData[["mzmlData"]]))
  
  temp2 <- mzmlData[["method"]][c(1, 2, 3, 4)]
  names(temp2) <- c("Identified metabolite", "Precursor m/z", "Product m/z", "Retention time")
  res <- merge(temp2, temp1, by.x = "Identified metabolite", by.y = "name", all.y = T)
  res[is.na(res)] <- 0
  
  res["sum"] <- rowSums(res[-c(1:4)])
  # colnames(res)
  res1 <- res %>%
    arrange(desc(sum))
  # len <- dim(res1)[1] * 0.95
  len <- dim(res1)[1]
  sum_col <- dim(res1)[2]
  res1 <- res1[1:len, -sum_col]
  res2 <- as.data.frame(t(res1))
  
  mzmlData[["peakRes"]] <- res2
  return(mzmlData)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "质谱数据",required = T)
  parser$add_argument("-o","--outputpath",default = "../raw/搜库数据", help = "保存路径")
  parser$add_argument("-m","--mode", help = "保存路径",required = T,choices = c("neg","pos"))
  parser$add_argument("-mp","--methodpath",default = "database/qualitative/QtargetLipid",
                      help = "保存路径")
  
  args <- parser$parse_args()
  
  args$methodpath <- selectfile(path = args$methodpath,file = paste0("qualitative-",args$mode,".txt"))
  
  result <- do.call(what = QtargetQualitative,args = args)
  
}
