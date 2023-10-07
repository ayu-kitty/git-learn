#' readregistration
#'
#' 读取内部分析单信息
#'
#' @param name 内部分析单文件名
#'
#' @export
readregistration <- function(name = "内部分析单.xlsx") {
  if (!file.exists(name)) {
    stop(paste0("不存在", name, "文件"))
  }

  data <- readdata(filename = name, sheet = "分析基本信息")
  data1 <- readdata(filename = name, sheet = "样本基本信息")
  data2 <- readdata(filename = name, sheet = "样本分组信息")
  if ("差异比较信息" %in% getsheetname(name)) {
    data3 <- readdata(filename = name, sheet = "差异比较信息")
  } else {
    data3 <- data.frame()
  }

  registration <- getregistration(data,
                                  data1,
                                  data2,
                                  data3)
  return(registration)
}

#' @export
getregistration <- function(data,
                            data1,
                            data2,
                            data3) {
  basic <- list(
    #2022.7.19新增客户名称和联系人名称读取-yxl
    "项目编号" = data[data[, 1] == "项目编号", 2],
    "客户单位" = data[data[, 1] == "客户单位", 2],
    "客户名称" = data[data[, 1] == "客户名称", 2],
    "联系人名称" = data[data[, 1] == "联系人", 2],
    "客户经理" = data[data[, 1] == "客户经理", 2],
    "客户邮箱" = data[data[, 1] == "客户邮箱", 2],
    "样本物种" = data[data[, 1] == "样本物种", 2],
    "样本类型" = data[data[, 1] == "样本类型", 2],
    "species" = NULLtoNA(data[data[, 1] == "映射物种", 2], nastop = T),
    "项目类型" = if (length(data[data[, 1] == "项目类别", 2]) == 0) {
      data[data[, 1] == "项目类型", 2]
    } else {
      data[data[, 1] == "项目类别", 2]
    },
    "处理类别" = paste0(
      data[data[, 1] == "QC质控", 2], "QC",
      data[data[, 1] == "比较分析", 2], "分析"
    ),
    "QC质控" = ifelse(data[data[, 1] == "QC质控", 2] == "有", T, F),
    "比较分析" = ifelse(data[data[, 1] == "比较分析", 2] == "有", T, F),
    "S3" = NULLtoNA(data[data[, 1] == "S3", 2], nastop = T),
    "分析编号" = NULLtoNA(data[data[, 1] == "分析编号", 2]),
    "项目报告" = NULLtoNA(data[data[, 1] == "项目报告", 2], fill = "项目报告"),
    "索引" = list(
      "ID" = NULLtoNA(data[data[, 1] == "ID", 2], nastop = T),
      "Metabolites" = NULLtoNA(data[data[, 1] == "Metabolites", 2], nastop = T),
      "Compound ID" = NULLtoNA(data[data[, 1] == "Compound ID", 2]),
      "Ion mode" = NULLtoNA(data[data[, 1] == "Ion mode", 2]),
      "KEGG" = NULLtoNA(data[data[, 1] == "KEGG", 2])
    )
  )

  handle <- list(
    missvalue = NULLtoNA(as.numeric(data[data[, 1] == "missvalue", 2])),
    rsd = NULLtoNA(as.numeric(data[data[, 1] == "rsd", 2])),
    zeroprocess = NULLtoNA(data[data[, 1] == "zeroprocess", 2]),
    samplenormalization = NULLtoNA(data[data[, 1] == "samplenormalization", 2]),
    datatransformation = NULLtoNA(data[data[, 1] == "datatransformation", 2]),
    datascaling = NULLtoNA(data[data[, 1] == "datascaling", 2]),
    log10L = NULLtoNA(as.logical(data[data[, 1] == "log10L", 2]), fill = F),
    PCAscaleC = NULLtoNA(data[data[, 1] == "PCAscaleC", 2], fill = "standard")
  )

  datafrom <- list(
    "LC原始文件" = list(
      negID = NULLtoNA(data[data[, 1] == "negID", 2]),
      posID = NULLtoNA(data[data[, 1] == "posID", 2]),
      negM = NULLtoNA(data[data[, 1] == "negM", 2]),
      posM = NULLtoNA(data[data[, 1] == "posM", 2])
    ),
    "GC原始文件" = NULLtoNA(data[data[, 1] == "GC原始文件", 2]),
    "其他原始文件" = NULLtoNA(data[data[, 1] == "数据矩阵", 2])
  )

  sample <- list(
    rawname = NULLtoNA(data1[, names(data1) == "实验名称"], nastop = F),
    samplename = NULLtoNA(data1[, names(data1) == "样本分析名称"], nastop = T),
    combat = NULLtoNA(data1[, names(data1) == "批次"]),
    deal = NULLtoNA(data1[, names(data1) == "处理"]),
    weight = NULLtoNA(data1[, names(data1) == "称重/体积"]),
    GCMS = list(
      raw = NULLtoNA(data1[, names(data1) == "GCMS-raw"]),
      mzml = NULLtoNA(data1[, names(data1) == "GCMS-mzml"]),
      name = NULLtoNA(data1[, names(data1) == "GCMS-name"])
    ),
    LCMS = list(
      neg = list(
        raw = NULLtoNA(data1[, names(data1) == "LCMS-neg-raw"]),
        mzml = NULLtoNA(data1[, names(data1) == "LCMS-neg-mzml"]),
        name = NULLtoNA(data1[, names(data1) == "LCMS-neg-name"])
      ),
      pos = list(
        raw = NULLtoNA(data1[, names(data1) == "LCMS-pos-raw"]),
        mzml = NULLtoNA(data1[, names(data1) == "LCMS-pos-mzml"]),
        name = NULLtoNA(data1[, names(data1) == "LCMS-pos-name"])
      )
    )
  )

  sample$class <- NULLtoNA(list(data1[, names(data1) == "分组"]), nastop = T)

  class <- list(
    group = NULLtoNA(data2[, names(data2) == "group"], nastop = T),
    samples = NULLtoNA(data2[, names(data2) == "samples"], nastop = T),
    type = NULLtoNA(data2[, names(data2) == "type"], nastop = T),
    shape = NULLtoNA(data2[, names(data2) == "shape"], nastop = T),
    fill = NULLtoNA(data2[, names(data2) == "fill"], nastop = T),
    colour = NULLtoNA(data2[, names(data2) == "colour"], nastop = T)
  )

  if (basic$`比较分析`) {
    if (nrow(data3[data3$class == "Group",]) == 0) {
      stop("请确认基本分组是否无比较分析")
    }

    dif <- list(
      compare = list(data3[data3$class == "Group", names(data3) == "compare"]),
      log10L = list(as.logical(data3[data3$class == "Group", names(data3) == "log10L"])),
      PCAscaleC = list(data3[data3$class == "Group", names(data3) == "PCAscaleC"]),
      PLSscaleC = list(data3[data3$class == "Group", names(data3) == "PLSscaleC"]),
      OPLSscaleC = list(data3[data3$class == "Group", names(data3) == "OPLSscaleC"]),
      PLSpermI = list(data3[data3$class == "Group", names(data3) == "PLSpermI"]),
      OPLSpermI = list(data3[data3$class == "Group", names(data3) == "OPLSpermI"]),
      UnivariateAnalysis = list(data3[data3$class == "Group", names(data3) == "UnivariateAnalysis"]),
      paired = list(data3[data3$class == "Group", names(data3) == "paired"]),
      p.adjust.method = list(data3[data3$class == "Group", names(data3) == "p.adjust.method"]),
      VIP = list(data3[data3$class == "Group", names(data3) == "VIP"]),
      FC = list(data3[data3$class == "Group", names(data3) == "FC"]),
      Pvalue = list(data3[data3$class == "Group", names(data3) == "Pvalue"]),
      adjPvalue = list(data3[data3$class == "Group", names(data3) == "adjPvalue"]),
      errorFC = list(data3[data3$class == "Group", names(data3) == "errorFC"]),
      order = list(data3[data3$class == "Group", names(data3) == "order"])
    )

    if (any(grepl(pattern = "^ExtendGroup[0-9]*", x = colnames(data1)))) {
      for (i in grep(pattern = "^ExtendGroup[0-9]*", x = colnames(data1))) {
        if (nrow(data3[data3$class == colnames(data1)[i], ]) == 0) {
          warning(paste0(colnames(data1)[i], "分组无比较分析"))
          next
        }
        sample$class <- c(sample$class, list(data1[, i]))
        dif$compare <- c(dif$compare, list(data3[data3$class == colnames(data1)[i], names(data3) == "compare"]))
        dif$log10L <- c(dif$log10L, list(as.logical(data3[data3$class == colnames(data1)[i], names(data3) == "log10L"])))
        dif$PCAscaleC <- c(dif$PCAscaleC, list(data3[data3$class == colnames(data1)[i], names(data3) == "PCAscaleC"]))
        dif$PLSscaleC <- c(dif$PLSscaleC, list(data3[data3$class == colnames(data1)[i], names(data3) == "PLSscaleC"]))
        dif$OPLSscaleC <- c(dif$OPLSscaleC, list(data3[data3$class == colnames(data1)[i], names(data3) == "OPLSscaleC"]))
        dif$PLSpermI <- c(dif$PLSpermI, list(data3[data3$class == colnames(data1)[i], names(data3) == "PLSpermI"]))
        dif$OPLSpermI <- c(dif$OPLSpermI, list(data3[data3$class == colnames(data1)[i], names(data3) == "OPLSpermI"]))
        dif$UnivariateAnalysis <- c(dif$UnivariateAnalysis, list(data3[data3$class == colnames(data1)[i], names(data3) == "UnivariateAnalysis"]))
        dif$paired <- c(dif$paired, list(data3[data3$class == colnames(data1)[i], names(data3) == "paired"]))
        dif$p.adjust.method <- c(dif$p.adjust.method, list(data3[data3$class == colnames(data1)[i], names(data3) == "p.adjust.method"]))
        dif$VIP <- c(dif$VIP, list(data3[data3$class == colnames(data1)[i], names(data3) == "VIP"]))
        dif$FC <- c(dif$FC, list(data3[data3$class == colnames(data1)[i], names(data3) == "FC"]))
        dif$Pvalue <- c(dif$Pvalue, list(data3[data3$class == colnames(data1)[i], names(data3) == "Pvalue"]))
        dif$adjPvalue <- c(dif$adjPvalue, list(data3[data3$class == colnames(data1)[i], names(data3) == "adjPvalue"]))
        dif$errorFC <- c(dif$errorFC, list(data3[data3$class == colnames(data1)[i], names(data3) == "errorFC"]))
        dif$order <- c(dif$order, list(data3[data3$class == colnames(data1)[i], names(data3) == "order"]))
      }
    }
  } else {
    dif <- NULL
  }

  registration <- list("info" = list("basic" = basic,
                                     "datafrom" = datafrom,
                                     "sample" = sample,
                                     "class" = class,
                                     "handle" = handle,
                                     "dif" = dif))

  attr(registration, "class") <- basic$S3
  return(registration)
}
