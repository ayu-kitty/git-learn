#!/opt/conda/bin/Rscript

#' 多元统计分析
#' 
#' @param data 数据
#' @param class 样本分组
#' @param group 分组组名
#' @param mode 多元统计模式
#' @param scaleC 标准化模式
#' @param predI 主成分数量
#' @param orthoI 正交主成分数量
#' @param permI 响应排序检验次数
#' @param log10L 是否log10处理
#' @param amount 最少主成分数量
#' @param adjust 对于响应排序检验结果进行调整参数运行
#' @param i 运行次数
#' @param maxi 最大运行次数
#'
#' @return list 
#' @export
mulstatistics <- function(data, 
                          class, 
                          group = "",
                          mode = "PCA",
                          scaleC = "standard",
                          predI = NA,
                          orthoI = 0,
                          permI = 0,
                          log10L = F,
                          amount = 3,
                          adjust = F,
                          i = 1,
                          maxi = 4) {
  
  options(warn = -1)
  
  if(any(is.na(data))){
    data[is.na(data)] <- min(as.matrix(data),na.rm = T)
  }

  if(is.na(permI)){permI <- 0}else if(permI < 0){permI <- 0}
  if(!is.na(predI)){if(predI < 0){predI <- NA}}
  if(!is.na(orthoI)){if(orthoI < 0){orthoI <- NA}}
  
  summarydata <- data.frame(Group = group,
                            Type = NA,
                            PRE = NA,
                            ORT = NA,
                            N = NA,
                            "R2X(cum)" = NA,
                            "R2Y(cum)" = NA,
                            "Q2(cum)" = NA,
                            R2 = NA,
                            Q2 = NA,
                            check.names = F,
                            stringsAsFactors = F)
  
  if (tolower(mode) == "pca") {
    if(!is.na(predI)){
      if(predI < amount){
        predI <- amount
      }
    }
    
    orthoI <- 0
    permI <- 0
    
    try <- try(
      {
        statistics <- ropls::opls(data,
                                  scaleC = scaleC,
                                  predI = predI,
                                  orthoI = orthoI,
                                  fig.pdfC = "none",
                                  info.txtC = "none",
                                  # printL = F,
                                  # plotL = F,
                                  log10L = log10L,
                                  algoC = "nipals",
                                  crossvalI = ifelse(dim(data)[1] < 7, dim(data)[1], 7)
        )
        model <- statistics@modelDF
      },
      silent = F
    )
    if (class(try) == "try-error") {
      statistics <- ropls::opls(data,
                                scaleC = scaleC,
                                predI = 2,
                                orthoI = orthoI,
                                fig.pdfC = "none",
                                info.txtC = "none",
                                # printL = F,
                                # plotL = F,
                                log10L = log10L,
                                algoC = "nipals",
                                crossvalI = ifelse(dim(data)[1] < 7, dim(data)[1], 7)
      )
      model <- statistics@modelDF
    }
    if(dim(statistics@summaryDF)[1] == 0){
      statistics@summaryDF[1, "pre"] <- 0
    }
    if (statistics@summaryDF[, "pre"] < amount) {
      statistics <- ropls::opls(data,
                                scaleC = scaleC,
                                predI = amount,
                                orthoI = orthoI,
                                fig.pdfC = "none",
                                info.txtC = "none",
                                # printL = F,
                                # plotL = F,
                                log10L = log10L,
                                algoC = "nipals",
                                crossvalI = ifelse(dim(data)[1] < 7, dim(data)[1], 7)
      )
      model <- statistics@modelDF
    }
  } else if (tolower(mode) == "pls" | tolower(mode) == "pls-da") {
    
    if(!is.na(predI)){
      if(predI < amount){
        predI <- amount
      }
    }
    
    orthoI <- 0
    
    try <- try(
      {
        statistics <- ropls::opls(
          x = data, y = class,
          scaleC = scaleC,
          predI = predI,
          orthoI = orthoI,
          permI = permI,
          fig.pdfC = "none",
          info.txtC = "none",
          # printL = F,
          # plotL = F,
          log10L = log10L,
          crossvalI = ifelse(dim(data)[1] < 7, dim(data)[1], 7)
        )
        model <- statistics@modelDF
      },
      silent = F
    )
    if (class(try) == "try-error") {
      statistics <- ropls::opls(
        x = data, y = class,
        scaleC = scaleC,
        predI = amount,
        orthoI = orthoI,
        permI = permI,
        fig.pdfC = "none",
        info.txtC = "none",
        # printL = F,
        # plotL = F,
        log10L = log10L,
        crossvalI = ifelse(dim(data)[1] < 7, dim(data)[1], 7)
      )
      model <- statistics@modelDF
    }
    if(dim(statistics@summaryDF)[1] == 0){
      statistics@summaryDF[1, "pre"] <- 0
    }
    if (statistics@summaryDF[, "pre"] < amount) {
      statistics <- ropls::opls(
        x = data, y = class,
        scaleC = scaleC,
        predI = amount,
        orthoI = orthoI,
        permI = permI,
        fig.pdfC = "none",
        info.txtC = "none",
        # printL = F, 
        # plotL = F,
        log10L = log10L,
        crossvalI = ifelse(dim(data)[1] < 7, dim(data)[1], 7)
      )
      model <- statistics@modelDF
    }
  } else if (tolower(mode) == "opls" | tolower(mode) == "opls-da") {
    
    if(!is.na(orthoI)){
      if(orthoI+1 < amount){
        orthoI <- amount-1
      }
    }
    
    predI <- 1
    
    if(length(unique(class)) > 2){
      warning("OPLS分析不支持3组及以上分析,现进行PLS-DA分析",immediate. = T)
      
      info <- mulstatistics(data = data, class = class, group = group, 
                            mode = "PLS-DA",
                            scaleC = scaleC, predI = NA, orthoI = 0, 
                            permI = permI, log10L = log10L,
                            amount = amount, adjust = adjust, i = 1, maxi = maxi)
      
      return(info)
    }
    
    try <- try(
      {
        statistics <- ropls::opls(
          x = data, y = class,
          scaleC = scaleC,
          predI = predI,
          orthoI = orthoI,
          permI = permI,
          fig.pdfC = "none",
          info.txtC = "none",
          # printL = F,
          # plotL = F,
          log10L = log10L,
          crossvalI = ifelse(dim(data)[1] < 7, dim(data)[1], 7)
        )
        model <- statistics@modelDF
      },
      silent = F
    )
    if (class(try) == "try-error") {
      statistics <- ropls::opls(
        x = data, y = class,
        scaleC = scaleC,
        predI = predI,
        orthoI = amount - 1,
        permI = permI,
        fig.pdfC = "none",
        info.txtC = "none",
        # printL = F,
        # plotL = F,
        log10L = log10L,
        crossvalI = ifelse(dim(data)[1] < 7, dim(data)[1], 7)
      )
      model <- statistics@modelDF
    }
    if(dim(statistics@summaryDF)[1] == 0){
      statistics@summaryDF[1, "pre"] <- 0
      statistics@summaryDF[1, "ort"] <- 0
    }
    if ((statistics@summaryDF[, "pre"] + statistics@summaryDF[, "ort"]) < amount) {
      statistics <- ropls::opls(
        x = data, y = class,
        scaleC = scaleC,
        predI = predI,
        orthoI = amount - 1,
        permI = permI,
        fig.pdfC = "none",
        info.txtC = "none",
        # printL = F,
        # plotL = F,
        log10L = log10L,
        crossvalI = ifelse(dim(data)[1] < 7, dim(data)[1], 7)
      )
      model <- statistics@modelDF
    }
  }
  
  if ((permI != 0) & (tolower(mode) == "opls" | tolower(mode) == "opls-da" | tolower(mode) == "pls" | tolower(mode) == "pls-da")) {
    print("响应排序检测开始")
    premdata <- as.data.frame(statistics@suppLs$permMN[, c(7, 2, 3)])
    names(premdata) <- c("ID", "R2", "Q2")
    premdata[, 1] <- abs(premdata[, 1] * 2 - 1)
    premdata1 <- data.frame(premdata[, c(1, 2)], class = "R2")
    names(premdata1)[2] <- c("va")
    maxpremdata11 <- premdata1[1, 1]
    maxpremdata12 <- premdata1[1, 2]
    premdata11 <- premdata1
    premdata11[, 1] <- premdata11[, 1] - maxpremdata11
    premdata11[, 2] <- premdata11[, 2] - maxpremdata12
    y <- premdata11[, 2]
    x <- premdata11[, 1]
    lmpre <- lm(y ~ x - 1)
    a1 <- round(lmpre$coefficients, 3)
    b1 <- round(maxpremdata12 - maxpremdata11 * a1, 3)
    
    premdata2 <- data.frame(premdata[, c(1, 3)], class = "Q2")
    names(premdata2)[2] <- c("va")
    maxpremdata21 <- premdata2[1, 1]
    maxpremdata22 <- premdata2[1, 2]
    premdata21 <- premdata2
    premdata21[, 1] <- premdata21[, 1] - maxpremdata21
    premdata21[, 2] <- premdata21[, 2] - maxpremdata22
    y <- premdata21[, 2]
    x <- premdata21[, 1]
    lmpre <- lm(y ~ x - 1)
    a2 <- round(lmpre$coefficients, 3)
    b2 <- round(maxpremdata22 - maxpremdata21 * a2, 3)
    premdata5 <- rbind(premdata1, premdata2)
    text2 <- paste0("R2=(0.0,", b1, ")   Q2=(0.0,", b2, ")")
    if (b2 > 0 | a2 < 0) {
      if (adjust & i <= maxi) {
        print("响应排序检测不符合要求，重新计算")
        
        if (a2 + b2 < 0) {
          x <- 0
          q2 <- a2 + b2
          q2_1 <- -1
          while (q2_1 <= q2 & x <= 4) {
            x <- x + 1
            q2_1 <- q2
            q2 <- ropls::opls(
              x = data, y = class,
              scaleC = scaleC,
              predI = predI,
              orthoI = x,
              permI = permI,
              fig.pdfC = "none",
              info.txtC = "none",
              # printL = F,
              # plotL = F,
              log10L = log10L,
              crossvalI = ifelse(dim(data)[1] < 7, dim(data)[1], 7)
            )@summaryDF[, "Q2(cum)"]
          }
          orthoI1 <- x
          info <- mulstatistics(
            data = data, class = class, group = group, mode = mode,
            scaleC = scaleC, predI = predI, orthoI = orthoI1, permI = permI, log10L = log10L,
            amount = amount, adjust = adjust, i = maxi + 1, maxi = maxi
          )
        } else {
          orthoI1 <- ort(orthoI = statistics@summaryDF[, "ort"], q2 = a2 + b2, b2 = b2)
          info <- mulstatistics(
            data = data, class = class, group = group, mode = mode,
            scaleC = scaleC, predI = predI, orthoI = orthoI1, permI = permI, log10L = log10L,
            amount = amount, adjust = adjust, i = i + 1, maxi = maxi
          )
        }
        return(info)
      }
      
      if (b2 > 0) print(paste0(group, "组响应排序检测Q2交集>0,请手动调整参数运算"))
      if (a2 < 0) print(paste0(group, "组响应排序检测Q2<Q2截距,请手动调整参数运算"))
    }
  }
  
  summary1 <- statistics@descriptionMC
  summary2 <- ropls::getSummaryDF(statistics)
  summarydata[, "Type"] <- mode
  summarydata[, "PRE"] <- summary2[, "pre"]
  summarydata[, "ORT"] <- summary2[, "ort"]
  summarydata[, "N"] <- as.numeric(summary1[1, 1])
  summarydata[, "R2X(cum)"] <- summary2[, "R2X(cum)"]
  if (tolower(mode) == "opls" | tolower(mode) == "opls-da" | tolower(mode) == "pls" | tolower(mode) == "pls-da") {
    summarydata[, "R2Y(cum)"] <- summary2[, "R2Y(cum)"]
    summarydata[, "Q2(cum)"] <- summary2[, "Q2(cum)"]
    if (permI != 0) {
      summarydata[, "R2"] <- b1
      summarydata[, "Q2"] <- b2
    }
  }
  
  if (tolower(mode) == "opls" | tolower(mode) == "opls-da") {
    y1 <- "to[1]"
  } else {
    y1 <- "t[2]"
  }
  
  if ((permI != 0) & (tolower(mode) == "opls" | tolower(mode) == "opls-da" | tolower(mode) == "pls" | tolower(mode) == "pls-da")) {
    info <- list(
      statistics = statistics,
      summarydata = summarydata,
      prem = list(data = premdata5, text = text2, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
    )
  } else {
    info <- list(
      statistics = statistics,
      summarydata = summarydata
    )
  }
  
  # options(warn=1)
  print(paste0("~",paste0(group,collapse = "-vs-"),"的",mode,"模式多元统计计算完成"))
  return(info)
}

#' @export
ort <- function(orthoI, q2, b2) {
  print(paste0("正交主成分数量：", orthoI))
  print(paste0("Q2：", q2))
  print(paste0("Q2截距：", b2))
  
  if (q2 < 0.6) {
    orthoI1 <- orthoI + 1
  } else if (q2 <= 0.9 & b2 < 0.1) {
    orthoI1 <- orthoI
  } else if (b2 > 0.4) {
    orthoI1 <- orthoI - 1
  } else if (q2 <= 0.9) {
    orthoI1 <- orthoI + 1
  } else if (q2 > 0.9) {
    orthoI1 <- orthoI - 1
  } else {
    orthoI1 <- orthoI
  }
  
  if (orthoI1 > 5) {
    orthoI1 <- 5
  }
  print(paste0("更新正交主成分数量：", orthoI1))
  return(orthoI1)
}

#' 运行多元统计分析
#' 
#' @param datafile 数据路径
#' @param classfile 分组路径
#' @param mulstatistics_rds_savepath 结果保存路径
#' @param ... 见`mulstatistics`
#'
#' @export
mulstatistics_file <- function(datafile = "oecloud/rawdata/datafile.txt",
                               classfile = "oecloud/rawdata/classfile.yaml",
                               mulstatistics_rds_savepath = NULL,
                               name = "",
                               group,
                               both = "F",
                               log10L = F,
                               ...){
  
  data <- readdata(filename = datafile,row.names = 1)
  class <- readdata(filename = classfile)
  
  if(both == "T"){
    lcdata <- data[grepl(pattern = "~LCMS$",x = row.names(data)),]
    lcmulstatistics_rds_savepath <- if(is.null(mulstatistics_rds_savepath)){mulstatistics_rds_savepath}else{gsub(pattern = ".rds$",replacement = "~LCMS.rds",x = mulstatistics_rds_savepath)}
    lcresult <- mulstatistics_file(datafile = lcdata,
                                   classfile = class,
                                   mulstatistics_rds_savepath = lcmulstatistics_rds_savepath,
                                   name = paste0(name,"~LCMS"),
                                   group = group,
                                   both = "F",
                                   log10L = F,
                                   ...)
    gcdata <- data[grepl(pattern = "~GCMS$",x = row.names(data)),]
    gcmulstatistics_rds_savepath <- if(is.null(mulstatistics_rds_savepath)){mulstatistics_rds_savepath}else{gsub(pattern = ".rds$",replacement = "~GCMS.rds",x = mulstatistics_rds_savepath)}
    gcresult <- mulstatistics_file(datafile = gcdata,
                                   classfile = class,
                                   mulstatistics_rds_savepath = gcmulstatistics_rds_savepath,
                                   name = paste0(name,"~GCMS"),
                                   group = group,
                                   both = "F",
                                   log10L = T,
                                   ...)
    result <- list(lcresult=lcresult,gcresult=gcresult)
    class(result) <- "mulstatistics-both"
  }else{
    # 数据处理
    anaclass <- class[group]
    if(any(duplicated(Reduce(c,anaclass)))){
      stop(paste0(name,"组中有重复样本"))
    }
    singleclass <- data.frame(Group = NULL)
    for ( i in 1:length(anaclass)) {
      newclass <- anaclass[[i]]
      if(!all(newclass %in% colnames(data))){
        warning(paste0("注意:",paste(newclass[!(newclass %in% colnames(data))],collapse = ";"),"不在数据中"),immediate. = T)
        newclass <- newclass[newclass %in% colnames(data)]
      }
      newsingleclass <- data.frame(Group = rep(names(anaclass)[i],length(newclass)),
                                   check.names = F,
                                   stringsAsFactors = F)
      rownames(newsingleclass) <- newclass
      singleclass <- rbind(singleclass,newsingleclass)
    }
    class <- singleclass[,1]
    data <- data[,row.names(singleclass),drop=F]
    data <- t(data)
    
    # 多元统计运算
    tryrun <- try({
      result <- mulstatistics(data = data, 
                              class = class,
                              group = name,
                              log10L = log10L,
                              ...)
    },silent = T)
    
    if("try-error" %in% class(tryrun)){
      
      result <- list(statistics = NULL,
                     summarydata = NULL)
      
    }
    
    result[["class"]] <- singleclass
    result[["group"]] <- factor(x = class,levels = unique(class))
    result[["groupname"]] <- rownames(data)
    class(result) <- "mulstatistics"
  }
  
  # print(paste0("~",paste0(group,collapse = "-vs-"),"多元统计计算完成"))
  
  if(is.null(mulstatistics_rds_savepath)){
    return(result)
  }else{
    # 数据保存
    saverds(data = result,
            filename = mulstatistics_rds_savepath)
    return(result)
  }
  
}

#' 运行多元统计分析，匹配meta包中流程
#' 
#' @param obj 数据
#' @param datafile 数据路径
#' @param classfile 分组路径
#' @param mulstatisticsyamlfile 多元统计yaml保存路径
#' @param mulstatistics_rds_savepath 结果保存路径
#' @param cores 核心数
#'
#' @export
run_mulstatisticsana <- function(obj,
                                 datafile = "oecloud/rawdata/datafile.txt",
                                 classfile = "oecloud/rawdata/classfile.yaml",
                                 mulstatisticsyamlfile = "oecloud/rawdata/mulstatistics.yaml",
                                 mulstatistics_rds_savepath = "oecloud/mulstatisticsanalyst",
                                 cores = 5,
                                 all = T){
  organizeobj(obj = obj,
              datafile = datafile,
              classfile = classfile)
  
  writemulstatisticsyaml(obj = obj,
                         datafile = datafile,
                         classfile = classfile,
                         mulstatisticsyamlfile = mulstatisticsyamlfile,
                         all = all)
  
  system(paste0("source /etc/profile;snakemake",
                " --snakefile ",packagepath(path = "snakemake/MultivariateStatisticalAnalysis-2023/mulstatisticana/whole/mulstatisticana_whole.smk"),
                " --directory ./",
                " --configfile ",mulstatisticsyamlfile,
                " --cores ",cores,
                " --quiet ",
                " --keep-going"),
         ignore.stdout = F, ignore.stderr = F)
  
  obj <- getmulstatistics_rdsforobj(obj = obj,
                                    mulstatistics_rds_savepath = mulstatistics_rds_savepath,
                                    all = all)
  
  unlink(x = "oecloud",recursive = T)
  unlink(x = ".snakemake",recursive = T)
  unlink(x = "classtype.xlsx")
  
  return(obj)
}
