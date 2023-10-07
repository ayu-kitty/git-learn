#!/opt/conda/bin/Rscript

#' M数据提取
#' 
#' @param code.negM 负离子M文件
#' @param code.posM 正离子M文件
#' @param negname 负离子M中样本名称 
#' @param posname 正离子M中样本名称 
#' @param samplename 样本分析名称
#' @param mode 处理模式，默认为"QI自带归一化模式"
#' @param deal mode为`raw`时,处理方式默认为"峰面积归一化"
#' @param neglab deal为`内标归一化`时,负离子内标名称
#' @param poslab deal为`内标归一化`时,正离子内标名称
#'
#' @export
dataM <- function(code.negM, 
                  code.posM, 
                  negname = NULL, 
                  posname = NULL, 
                  samplename = NULL,
                  weight = NULL,
                  mode = "normal", 
                  deal = "峰面积归一化", 
                  neglab = NULL, 
                  poslab = NULL) {
  print("定量数据处理")
  print("neg定量数据提取")
  negM <- Mfetch(name = code.negM,
                 grouping = negname,
                 weight = weight,
                 samplename = samplename,
                 mode = mode, 
                 deal = deal, 
                 lab = neglab)
  
  print("pos定量数据提取")
  posM <- Mfetch(name = code.posM,
                 grouping = posname,
                 samplename = samplename,
                 weight = weight,
                 mode = mode, 
                 deal = deal, 
                 lab = poslab)
  
  print("正负离子数据合并")
  negposM <- negposMmerge(datanegM = negM,
                          dataposM = posM)
  
  print("定量数据处理完毕")
  return(negposM)
}

#' 获取M数据
#' 
#' @param name M文件名
#' @param grouping M中样本名称 
#' @param samplename 样本分析名称
#' @param mode 处理模式，normal默认为"QI自带归一化模式",raw为提取未归一化数据
#' @param deal mode为`raw`时,处理方式默认为"峰面积归一化","内标归一化","细胞数归一化"
#' @param lab deal为`内标归一化`时,内标名称
#'
#' @export
Mfetch <- function(name, 
                   grouping = NULL, 
                   samplename = NULL,
                   weight = NULL,
                   mode = "normal", 
                   deal = "峰面积归一化", 
                   lab = NULL) {
  tryfetch <- try(
    {
      data <- read.csv(name, header = TRUE, skip = 2, check.names = F, encoding = "UTF-8", stringsAsFactors = F)
      data2 <- read.csv(name, header = TRUE, check.names = F, nrows = 2)
      
      srow <- which(colnames(data2) == "Normalised abundance")
      rowm <- which(colnames(data2) == "Raw abundance") - 1
      data3 <- subset(data, select = c(Compound, `m/z`, `Retention time (min)`, Charge))
      
      if (mode == "normal") {
        data4 <- data[, srow:rowm, drop = F]
      } else if (mode == "raw") {
        data4 <- data[, (rowm + 1):(rowm - srow + rowm + 1), drop = F]
        if (deal == "峰面积归一化") {
          print("mode = raw,deal = 峰面积归一化")
          i <- 1
          while (i <= dim(data4)[2]) {
            data4[, i] <- data4[, i] / sum(data4[, i]) * 10000
            i <- i + 1
          }
        } else if (deal == "内标归一化") {
          print("mode = raw,deal = 内标归一化")
          if (!is.null(lab)) {
            lab2 <- data3[data3[, "Compound"] %in% lab, c("Compound", "Retention time (min)")]
            lab1 <- lab2[order(lab2[, 2]), 2]
            lab1 <- c(0, lab1, Inf)
            data5 <- data4
            i <- 1
            while (i < length(lab1)) {
              if (lab1[i + 1] == Inf) {
                num <- lab1[i]
              } else {
                num <- lab1[i + 1]
              }
              num <- lab2[lab2[, 2] == num, 1]
              
              j <- 1
              while (j <= dim(data4)[2]) {
                data4[lab1[i] < data3[, "Retention time (min)"] & data3[, "Retention time (min)"] <= lab1[i + 1], j] <- data4[lab1[i] < data3[, "Retention time (min)"] & data3[, "Retention time (min)"] <= lab1[i + 1], j] / data5[data3[, "Compound"] == num, j]
                j <- j + 1
              }
              i <- i + 1
            }
          }
        }else if (deal == "细胞数归一化") {
          print("mode = raw,deal = 细胞数归一化")
          if(any(is.na(as.numeric(weight))) == T){
            print("称重数据中含有非数字字符,例如单位,请检查")
          }else{
            dataweight <- t(data.frame(weight,row.names = samplename))
            dataweight <- dataweight[,colnames(data4)]
            data5 <- t(data4)/as.numeric(dataweight)
            data4 <- t(data5)
            data4 <- as.data.frame(data4)
          }
        }
      }
      
      if(!is.null(samplename)){
        if (all(grouping %in% colnames(data4))) {
          data4 <- data4[, grouping]
          colnames(data4) <- samplename
        } else {
          warning("样本量与样品信息不符合")
        } 
      }
      
      data <- cbind(data3, data4)
      colnames(data)[1] <- "ID"
      colnames(data)[4] <- "Ion mode"
    },
    silent = F
  )
  
  if (class(tryfetch) == "try-error") {
    data <- NULL
  }
  
  return(data)
}

#' negM和posM合并
#'
#' @param datanegM 负离子M数据
#' @param dataposM 正离子M数据
#'
#' @export
negposMmerge <- function(datanegM, 
                         dataposM) {
  negM <- datanegM
  posM <- dataposM
  if (!is.null(negM)) {
    negM$`Ion mode` <- "neg"
  }
  if (!is.null(posM)) {
    posM$`Ion mode` <- "pos"
  }
  
  # 将合并-neg-M.csv和-pos-M.csv数据
  negposM <- rbind(negM, posM)
  negposM <- negposM[!duplicated(negposM$ID), ]
  
  return(negposM)
}
