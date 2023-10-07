#' ArrangeInfoFrom
#'
#' 整理数据矩阵表格为项目登记单或者info信息进行后续分析
#'
#' @param name 默认为数据矩阵.xlsx
#' @param datafrom  默认为数据矩阵 sheet
#' @param indexfrom 默认为索引 sheet
#' @param groupfrom 默认为分组 sheet
#' @param comparefrom 默认为比较 sheet
#' @param ... 见[ArrangeInfoFromData()]
#'
#' @return 登记单列表
#'
#' @export
ArrangeInfoFrom <- function(name = "数据矩阵.xlsx",
                            datafrom = "数据矩阵",
                            indexfrom = "索引",
                            groupfrom = "分组",
                            comparefrom = "比较",
                            overwrite = T,
                            ...) {
  print("分析确认单整理中")

  sheetname <- openxlsx::getSheetNames(name)

  if ((datafrom %in% sheetname) &
    (indexfrom %in% sheetname) &
    (indexfrom %in% sheetname)) {
    data <- readdata(name = name, sheet = datafrom)
    index <- readdata(name = name, sheet = indexfrom)
    group <- readdata(name = name, sheet = groupfrom)
    if (comparefrom %in% sheetname) {
      compare <- readdata(name = name, sheet = comparefrom)
      
      if(!("分组" %in% colnames(compare))){
        compare[,"分组"] <- colnames(group)[2]
      }
      
    } else {
      warning(paste0(name, "文件中无", comparefrom, "-sheet"), immediate. = T)
      compare <- NULL
    }
  } else if (datafrom %in% sheetname) {
    if (indexfrom %in% sheetname) {
    } else {
      warning(paste0(name, "文件中无", indexfrom, "-sheet"), immediate. = T)
    }

    if (groupfrom %in% sheetname) {
    } else {
      warning(paste0(name, "文件中无", groupfrom, "-sheet"), immediate. = T)
    }
    print("使用其他模式处理数据")

    data <- readdata(name = name, sheet = datafrom)
    data1 <- data[, !is.na(data[1, ])]
    if (any(data1[1, ] == "ID") &
      any(data1[1, ] == "Metabolites")) {
      class1 <- as.vector(t(data[1, !is.na(data[1, ]) & data[1, ] != "ID" &
        data[1, ] != "Metabolites" &
        data[1, ] != "Compound ID" &
        data[1, ] != "Ion mode" &
        data[1, ] != "KEGG"]))
      rawname <- names(data1)[data1[1, ] %in% class1]
      index <- data1[1, !(data1[1, ] %in% class1)]
      index <- data.frame(
        "索引" = as.vector(t(index)[, 1]),
        "列名" = names(index),
        stringsAsFactors = F
      )
      group <- data.frame(
        sample = rawname,
        group = class1,
        stringsAsFactors = F
      )
      data <- data[-1, ]
      data[, colnames(data) %in% rawname] <- apply(data[, colnames(data) %in% rawname], 2, as.numeric)
    } else {
      stop("无ID或Metabolites的关键索引")
    }

    if (comparefrom %in% sheetname) {
      compare <- readdata(name = name, sheet = comparefrom)

      if("配对" %in% colnames(compare)){
        compare[is.na(compare[, "配对"]), "配对"] <- F
      }else{
        compare[, "配对"] <- F
      }

      compare[, "分组"] <- "group"
      print(compare)
    } else {
      warning(paste0(name, "文件中无", comparefrom, "-sheet"), immediate. = T)
      compare <- NULL
    }
  } else {
    if (datafrom %in% sheetname) {
    } else {
      stop(paste0(name, "文件中无", datafrom, "-sheet"))
    }
  }

  registration <- ArrangeInfoFromData(
    data = data,
    index = index,
    group = group,
    compare = compare,
    overwrite = overwrite,
    ...
  )

  return(registration)
}

#' ArrangeInfoFromData
#'
#' 获取相应相信返回info列表
#'
#' @param data 数据矩阵
#' @param index 索引
#' @param group 分组
#' @param compare 比较
#' @param type 类型，默认为非靶向代谢-外来数据
#' @param species 物种，默认为ko
#' @param missvalue 缺失值筛选，默认为0.5
#' @param rsd rsd筛选，默认为0.3
#' @param zeroprocess 0值填充，默认为min
#' @param samplenormalization 标准化
#' @param datatransformation log转化
#' @param datascaling 归一化
#' @param log10L 多元统计log10处理
#' @param PCAscaleC 多元统计PCA归一化模式
#' @param PLSscaleC 多元统计PLS归一化模式
#' @param OPLSscaleC 多元统计OPLS归一化模式
#' @param PLSpermI 多元统计PLS检验次数
#' @param OPLSpermI 多元统计OPLS检验次数
#' @param UnivariateAnalysis P值算法
#' @param p.adjust.method 矫正P值算法
#' @param VIP 多元统计VIP筛选标准，默认为1
#' @param FC FC筛选标准
#' @param Pvalue P值筛选标准，默认为0.05
#' @param adjPvalue 矫正P值筛选标准
#' @param errorFC 无P或VIP筛选情况下，启用FC筛选标准
#' @param order 排序方式，默认为VIP
#' @param saveregistration 逻辑，是否包含内部分析单.xlsx
#'
#' @return 登记单列表
#'
#' @export
ArrangeInfoFromData <- function(data,
                                index,
                                group,
                                compare,
                                type = "非靶向代谢-外来数据",
                                species = "ko",
                                missvalue = 0.5,
                                rsd = 0.3,
                                zeroprocess = "min",
                                samplenormalization = NA,
                                datatransformation = NA,
                                datascaling = NA,
                                log10L = FALSE,
                                PCAscaleC = "standard",
                                PLSscaleC = "pareto",
                                OPLSscaleC = "pareto",
                                PLSpermI = 0,
                                OPLSpermI = 200,
                                UnivariateAnalysis = "ttest",
                                p.adjust.method = "BH",
                                VIP = 1,
                                FC = NA,
                                Pvalue = 0.05,
                                adjPvalue = NA,
                                errorFC = 2,
                                order = "VIP",
                                saveregistration = T,
                                overwrite = T) {
  library(openxlsx)
  wb <- createWorkbook()

  addWorksheet(wb, "数据矩阵")
  writeData(wb, sheet = "数据矩阵", data)

  index <- index[!is.na(index[, 2]), ]
  index <- index[index[, 2] != "", ]
  index <- index[index[, 1] %in% c("ID", "Metabolites", "Compound ID", "Ion mode", "KEGG"), ]
  index <- index[, 1:2, drop = F]
  colnames(index) <- c("key", "value")

  if (any(index[, 1] == "ID") & any(index[, 1] == "Metabolites")) {
    if ((index[index[, 1] == "ID", 2] %in% colnames(data)) &
      (index[index[, 1] == "Metabolites", 2] %in% colnames(data))) {
      print("~~分析基本信息整理中")

      info <- data.frame(key = "项目报告", value = "项目报告", stringsAsFactors = F)
      info["项目类别", ] <- c("项目类别", type)
      info["映射物种", ] <- c("映射物种", species)
      info["S3", ] <- c("S3", "OutData")
      info["索引", ] <- c("索引", "")
      info <- rbind(info, index)

      info["预处理", ] <- c("预处理", "")
      info["missvalue", ] <- c("missvalue", missvalue)
      info["rsd", ] <- c("rsd", rsd)
      info["zeroprocess", ] <- c("zeroprocess", zeroprocess)
      info["samplenormalization", ] <- c("samplenormalization", samplenormalization)
      info["datatransformation", ] <- c("datatransformation", datatransformation)
      info["datascaling", ] <- c("datascaling", datascaling)
      info["分析参数", ] <- c("分析参数", "")
      info["log10L", ] <- c("log10L", log10L)
      info["PCAscaleC", ] <- c("PCAscaleC", PCAscaleC)

      info["其他原始文件", ] <- c("其他原始文件", "")
      info["数据矩阵", ] <- c("数据矩阵", "内部分析单.xlsx")

      group <- group[!is.na(group[, 1]), ]
      group <- group[!is.na(group[, 2]), ]
      sample <- group[, c(1, 2)]
      colnames(sample) <- c("样本分析名称", "分组")

      if (!all(sample[, 1] %in% colnames(data))) {
        stop("分组表中有分组信息的样本未找到对应数据列")
      }

      print("~~样本信息及分组信息整理中")
      sample[, "样本下机名称"] <- sample[, "样本分析名称"]
      sample[, "批次"] <- 1
      sample[, "处理"] <- NA

      samplegroup <- data.frame(group = unique(sample[, "分组"]), stringsAsFactors = F)
      for (i in 1:nrow(samplegroup)) {
        samplegroup[i, "samples"] <- paste(sample[sample[, "分组"] == samplegroup[i, "group"], "样本分析名称"],
          sep = "", collapse = ","
        )
        samplegroup[i, "type"] <- "常规"
      }

      if ("QC" %in% samplegroup[, 1]) {
        info["QC", ] <- c("QC质控", "有")
        samplegroup <- rbind(
          samplegroup[samplegroup[, "group"] == "QC", ],
          samplegroup[samplegroup[, "group"] != "QC", ]
        )
      } else {
        warning("~~~~未发现质控样本信息", immediate. = T)
        info["QC", ] <- c("QC质控", "无")
        info["rsd", ] <- c("rsd", NA)
      }

      if (!is.null(compare)) {
        print("~~差异比较信息整理中")
        if (all(colnames(compare) %in% c("比较组", "配对", "分组"))) {
          groupname <- colnames(group)[2]
          compare2 <- compare[compare[, "分组"] == groupname, ]
          if (nrow(compare2) == 0) {
            stop(paste0("差异比较信息中需提供", groupname, "的比较组"))
          } else {
            comparegroup <- data.frame(
              compare = compare2[, "比较组"],
              paired = compare2[, "配对"],
              class = "Group",
              stringsAsFactors = F
            )
            for (j in 1:nrow(comparegroup)) {
              classname <- unlist(strsplit(comparegroup[j, "compare"], split = "/"))
              if (all(classname %in% group[, 2])) {
                comparegroup[j, "num"] <- length(classname)
              } else {
                stop(paste0(groupname, "分组中未找到", comparegroup[j, "compare"], "相关样本"))
              }
            }
          }

          if (ncol(group) > 2) {
            for (i in 3:ncol(group)) {
              groupname <- colnames(group)[i]
              compare2 <- compare[compare[, "分组"] == groupname, ]
              if (nrow(compare2) == 0) {
                warning(paste0("差异比较信息中需提供", groupname, "的比较组"), immediate. = T)
                next
              } else {
                comparegroup2 <- data.frame(
                  compare = compare2[, "比较组"],
                  paired = compare2[, "配对"],
                  class = paste0("ExtendGroup", i - 2),
                  stringsAsFactors = F
                )
                for (j in 1:nrow(comparegroup2)) {
                  classname <- unlist(strsplit(comparegroup2[j, "compare"], split = "/"))
                  if (all(classname %in% group[, i])) {
                    comparegroup2[j, "num"] <- length(classname)
                  } else {
                    warning(paste0(groupname, "分组中未找到", comparegroup2[j, "compare"], "相关样本"))
                    next
                  }
                }

                sample[, paste0("ExtendGroup", i - 2)] <- group[, i]
                samplegroup2 <- data.frame(group = unique(unlist(strsplit(comparegroup2[, "compare"], split = "/"))), stringsAsFactors = F)
                for (j in 1:nrow(samplegroup2)) {
                  samplegroup2[j, "samples"] <- paste(sample[sample[, paste0("ExtendGroup", i - 2)] == samplegroup2[j, "group"], "样本分析名称"],
                    sep = "", collapse = ","
                  )
                  samplegroup2[j, "type"] <- "扩展"
                }
                samplegroup <- rbind(samplegroup, samplegroup2)
                if (any(duplicated(samplegroup[, "group"]))) {
                  warning(paste0(groupname, "分组与其他分组有相同组名"),immediate. = T)
                }
                comparegroup <- rbind(comparegroup, comparegroup2)
              }
            }
          }

          info["比较分析", ] <- c("比较分析", "有")
        } else {
          stop("差异比较信息中需提供比较组、配对、分组信息")
        }

        comparegroup[, "log10L"] <- log10L
        comparegroup[, "PCAscaleC"] <- PCAscaleC
        comparegroup[, "PLSscaleC"] <- PLSscaleC
        comparegroup[, "OPLSscaleC"] <- OPLSscaleC
        comparegroup[, "PLSpermI"] <- PLSpermI
        comparegroup[, "OPLSpermI"] <- OPLSpermI
        comparegroup[comparegroup$num > 2, "PLSpermI"] <- 200
        comparegroup[comparegroup$num > 2, "OPLSpermI"] <- 0
        comparegroup[, "UnivariateAnalysis"] <- UnivariateAnalysis
        comparegroup[, "p.adjust.method"] <- p.adjust.method

        comparegroup[, "VIP"] <- VIP
        comparegroup[, "FC"] <- FC
        comparegroup[, "Pvalue"] <- Pvalue
        comparegroup[, "adjPvalue"] <- adjPvalue
        comparegroup[, "errorFC"] <- errorFC
        comparegroup[, "order"] <- order

        addWorksheet(wb, "差异比较信息")
        writeData(wb, sheet = "差异比较信息", comparegroup)
      } else {
        comparegroup <- data.frame()
        warning("~~~~未发现差异比较信息", immediate. = T)
        info["比较分析", ] <- c("比较分析", "无")
      }

      if (dim(samplegroup)[1] <= 300) {
        samplegroup[, "shape"] <- rep(c(21, 22, 24, 23, 25), 60)[1:dim(samplegroup)[1]]
        # shape<-c("圆","方","上三","菱","下三")
        fill <- c(
          "green4", "blue3", "firebrick",
          "gold", "darkviolet", "darkorange",
          "skyblue3", "olivedrab3", "dodgerblue3",
          "aquamarine2", "deeppink3", "slateblue3",
          "brown2", "palegreen2", "chocolate2",
          "antiquewhite3", "steelblue1", "violetred1",
          "burlywood3", "pink1", "slategray2",
          "orangered1", "cyan3", "yellow4",
          "red", "plum", "greenyellow",
          "mediumpurple2", "tan1", "magenta"
        )
        if (dim(samplegroup)[1] <= 30) {
          fill <- fill[1:dim(samplegroup)[1]]
        } else {
          set.seed(111)
          colorrange <- colors()[c(1:150, 361:650)]
          colorrange <- colorrange[!(colorrange %in% fill)]
          fill <- c(fill, sample(colorrange, dim(samplegroup)[1] - length(fill)))
        }

        samplegroup[, "fill"] <- fill
        samplegroup[, "colour"] <- rep("dimgrey", 300)[1:dim(samplegroup)[1]]
      } else {
        stop("样本分组最多支持300组，现超过限制，请减少样本分组")
      }

      addWorksheet(wb, "分析基本信息")
      writeData(wb, sheet = "分析基本信息", info)
      addWorksheet(wb, "样本基本信息")
      writeData(wb, sheet = "样本基本信息", sample)
      addWorksheet(wb, "样本分组信息")
      writeData(wb, sheet = "样本分组信息", samplegroup)

      if (saveregistration) {
        saveWorkbook(wb, file = "内部分析单.xlsx", overwrite = overwrite)
      }
    } else {
      stop("数据中无ID或Metabolites的关键索引")
    }
  } else {
    stop("无ID或Metabolites的关键索引")
  }

  registration <- getregistration(
    data = info,
    data1 = sample,
    data2 = samplegroup,
    data3 = comparegroup
  )
  registration$data$rawdata$data <- data
  registration <- autodataprocess(registration, fetch = F)
  return(registration)
  print("complete!")
}
