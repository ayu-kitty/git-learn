#' @export
automissvalue <- function(data,
                          group = T,
                          missvalue = data$info$handle$missvalue) {
  data3 <- data
  data3$info$handle$missvalue <- missvalue
  
  if (!is.na(data3$info$handle$missvalue)) {
    print("缺失值处理开始")
    data1 <- data3$data$data$data
    
    if (is.na(data3$info$basic$索引$`Ion mode`)) {
      rawnum <- c(all = dim(data1)[1])
    } else {
      a <- data3$data$predata$information[, c(
        data3$info$basic$索引$ID,
        data3$info$basic$索引$`Ion mode`
      )]
      b <- tolower(whto(a, rownames(data1)))
      rawnum <- c(
        pos = length(b[b == "pos"]),
        neg = length(b[b == "neg"]),
        all = dim(data1)[1]
      )
    }
    
    data3$data$missvalue$num$rawnum <- rawnum
    
    if (group == T) {
      group1 <- unique(whto(
        data.frame(
          samplename = data3$info$sample$samplename,
          class = data3$info$sample$class[[1]],
          stringsAsFactors = F
        ),
        colnames(data1)
      ))
      data2 <- NA
      
      i <- 1
      while (i <= length(group1)) {
        if (group1[i] != "QC") {
          data4 <- data.frame(`default-value` = apply(data1[, (whto(
            data.frame(
              samplename = data3$info$sample$samplename,
              class = data3$info$sample$class[[1]],
              stringsAsFactors = F
            ),
            colnames(data1)
          ) %in% group1[i]), drop = F],
          MARGIN = 1, FUN = defaultvalue
          ), check.names = F)
          data2 <- cbind(data2, data4)
        }
        i <- i + 1
      }
      
      data2 <- data.frame(`default-value` = apply(data2[, -1, drop = F], MARGIN = 1, FUN = max), check.names = F)
    } else {
      data2 <- data.frame(`default-value` = apply(data1[, !(whto(
        data.frame(
          samplename = data3$info$sample$samplename,
          class = data3$info$sample$class[[1]], stringsAsFactors = F
        ),
        colnames(data1)
      ) %in% "QC"), drop = F],
      MARGIN = 1, FUN = defaultvalue
      ), check.names = F)
    }
    
    data3$data$missvalue$missvalue <- data2
    
    data1 <- data1[data2[, "default-value"] >= data3$info$handle$missvalue, ]
    data3$data$missvalue$data <- data3$data$data$data
    data3$data$missvalue$dealdata <- data1
    data3$data$data$data <- data1
    
    if (is.na(data3$info$basic$索引$`Ion mode`)) {
      num <- c(all = dim(data1)[1])
    } else {
      a <- data3$data$predata$information[, c(
        data3$info$basic$索引$ID,
        data3$info$basic$索引$`Ion mode`
      )]
      b <- tolower(whto(a, rownames(data1)))
      num <- c(
        pos = length(b[b == "pos"]),
        neg = length(b[b == "neg"]),
        all = dim(data1)[1]
      )
    }
    
    yield <- num / rawnum
    data3$data$missvalue$num$num <- num
    data3$data$missvalue$num$yield <- yield
    
    print(data.frame(
      "缺省值筛选" = names(data3$data$missvalue$num$rawnum),
      "筛选前" = data3$data$missvalue$num$rawnum,
      "筛选后" = data3$data$missvalue$num$num,
      "得率" = data3$data$missvalue$num$yield
    ))
    
    print("缺失值处理完毕")
  } else {
    warning("未进行缺失值处理，请确认", immediate. = T)
  }
  
  return(data3)
}


# 数据矩阵缺省值计算
#' @export
defaultvalue <- function(x) {
  k <- 1 - sum(x == 0) / length(x)
  return(k)
}

#' @export
autorsdfilter <- function(data,
                          rsd = data$info$handle$rsd) {
  data3 <- data
  data3$info$handle$rsd <- rsd
  
  if (!is.na(data3$info$handle$rsd)) {
    print("rsd筛选处理开始")
    data1 <- data3$data$data$data
    
    if (is.na(data3$info$basic$索引$`Ion mode`)) {
      rawnum <- c(all = dim(data1)[1])
    } else {
      a <- data3$data$predata$information[, c(
        data3$info$basic$索引$ID,
        data3$info$basic$索引$`Ion mode`
      )]
      b <- tolower(whto(a, rownames(data1)))
      rawnum <- c(
        pos = length(b[b == "pos"]),
        neg = length(b[b == "neg"]),
        all = dim(data1)[1]
      )
    }
    
    data3$data$rsdfilter$num$rawnum <- rawnum
    
    data2 <- data.frame(
      `mean-qc` = apply(data1[, (whto(
        data.frame(
          samplename = data3$info$sample$samplename,
          class = data3$info$sample$class[[1]], stringsAsFactors = F
        ),
        colnames(data1)
      ) %in% "QC")], MARGIN = 1, FUN = mean),
      `sd-qc` = apply(data1[, (whto(
        data.frame(
          samplename = data3$info$sample$samplename,
          class = data3$info$sample$class[[1]], stringsAsFactors = F
        ),
        colnames(data1)
      ) %in% "QC")], MARGIN = 1, FUN = sd),
      check.names = F
    )
    data2[, "rsd-qc"] <- data2$`sd-qc` / data2$`mean-qc`
    
    data3$data$rsdfilter$rsd <- data2
    
    data2 <- data2[!is.na(data2$`rsd-qc`), ]
    data2 <- data2[data2$`rsd-qc` <= data3$info$handle$rsd, ]
    data1 <- data1[row.names(data2), ]
    
    data3$data$rsdfilter$data <- data3$data$data$data
    data3$data$rsdfilter$dealdata <- data1
    data3$data$data$data <- data1
    
    if (is.na(data3$info$basic$索引$`Ion mode`)) {
      num <- c(all = dim(data1)[1])
    } else {
      a <- data3$data$predata$information[, c(
        data3$info$basic$索引$ID,
        data3$info$basic$索引$`Ion mode`
      )]
      b <- tolower(whto(a, rownames(data1)))
      num <- c(
        pos = length(b[b == "pos"]),
        neg = length(b[b == "neg"]),
        all = dim(data1)[1]
      )
    }
    
    yield <- num / rawnum
    data3$data$rsdfilter$num$num <- num
    data3$data$rsdfilter$num$yield <- yield
    
    print(data.frame(
      "RSD筛选" = names(data3$data$rsdfilter$num$rawnum),
      "筛选前" = data3$data$rsdfilter$num$rawnum,
      "筛选后" = data3$data$rsdfilter$num$num,
      "得率" = data3$data$rsdfilter$num$yield
    ))
    
    print("rsd筛选处理完毕")
  } else {
    warning("未进行rsd筛选处理，请确认", immediate. = T)
  }
  
  return(data3)
}
