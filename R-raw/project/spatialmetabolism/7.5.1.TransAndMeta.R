#' @export
getparameter <- function(transspot, metaspot, xresolution, yresolution) {
  if ((dim(transspot)[1] == dim(metaspot)[1]) & dim(transspot)[1] >= 3) {
    transspot2 <- transspot
    transspot2$x <- (transspot$x / 2 + 0.5) * 100
    transspot2$y <- transspot$y * sqrt(7500)

    metaspot2 <- metaspot
    metaspot2$x <- metaspot$x * xresolution
    metaspot2$y <- metaspot$y * xresolution

    combination <- combn(x = dim(transspot)[1], m = 2)

    trans <- F
    rawdata <- data.frame(transdistance = NA, metadistance = NA, ratio = NA, transangle = NA, metaangle = NA, angle = NA, A = combination[1, ], B = combination[2, ])
    for (i in 1:ncol(combination)) {
      rawdata[i, "transdistance"] <- sqrt((transspot2[combination[1, i], 1] - transspot2[combination[2, i], 1])^2 + (transspot2[combination[1, i], 2] - transspot2[combination[2, i], 2])^2)
      rawdata[i, "metadistance"] <- sqrt((metaspot2[combination[1, i], 1] - metaspot2[combination[2, i], 1])^2 + (metaspot2[combination[1, i], 2] - metaspot2[combination[2, i], 2])^2)
      rawdata[i, "ratio"] <- rawdata[i, "metadistance"] / rawdata[i, "transdistance"]
      rawdata[i, "transangle"] <- atan2((transspot2[combination[1, i], 2] - transspot2[combination[2, i], 2]), (transspot2[combination[1, i], 1] - transspot2[combination[2, i], 1])) / pi * 180
      rawdata[i, "metaangle"] <- atan2((metaspot2[combination[1, i], 2] - metaspot2[combination[2, i], 2]), (metaspot2[combination[1, i], 1] - metaspot2[combination[2, i], 1])) / pi * 180
      angle <- rawdata[i, "metaangle"] - rawdata[i, "transangle"]
      rawdata[i, "angle"] <- ifelse(angle < 0, angle + 360, angle)
    }

    if (max(rawdata[, "angle"]) - min(rawdata[, "angle"]) > 90) {
      trans <- T
      metaspot2$y <- -metaspot$y * xresolution
      rawdata <- data.frame(transdistance = NA, metadistance = NA, ratio = NA, transangle = NA, metaangle = NA, angle = NA, A = combination[1, ], B = combination[2, ])
      for (i in 1:ncol(combination)) {
        rawdata[i, "transdistance"] <- sqrt((transspot2[combination[1, i], 1] - transspot2[combination[2, i], 1])^2 + (transspot2[combination[1, i], 2] - transspot2[combination[2, i], 2])^2)
        rawdata[i, "metadistance"] <- sqrt((metaspot2[combination[1, i], 1] - metaspot2[combination[2, i], 1])^2 + (metaspot2[combination[1, i], 2] - metaspot2[combination[2, i], 2])^2)
        rawdata[i, "ratio"] <- rawdata[i, "metadistance"] / rawdata[i, "transdistance"]
        rawdata[i, "transangle"] <- atan2((transspot2[combination[1, i], 2] - transspot2[combination[2, i], 2]), (transspot2[combination[1, i], 1] - transspot2[combination[2, i], 1])) / pi * 180
        rawdata[i, "metaangle"] <- atan2((metaspot2[combination[1, i], 2] - metaspot2[combination[2, i], 2]), (metaspot2[combination[1, i], 1] - metaspot2[combination[2, i], 1])) / pi * 180
        angle <- rawdata[i, "metaangle"] - rawdata[i, "transangle"]
        rawdata[i, "angle"] <- ifelse(angle < 0, angle + 360, angle)
      }
    }
  } else {
    stop("空转与空代数量不对应，或少于三个")
  }

  meanratio <- mean(rawdata[, "ratio"])
  maxratio <- max(rawdata[, "ratio"])
  minratio <- min(rawdata[, "ratio"])
  sdratio <- sd(rawdata[, "ratio"])
  rsdratio <- sdratio / meanratio

  if (maxratio > 1.05 | minratio < 0.95 | rsdratio > 0.1) {
    warning("缩放比不符合标准，请查询原因", immediate. = T)
  }

  meanangle <- mean(rawdata[, "angle"])
  maxangle <- max(rawdata[, "angle"])
  minangle <- min(rawdata[, "angle"])
  sdangle <- sd(rawdata[, "angle"])
  rsdangle <- sdangle / meanangle

  if (rsdangle > 0.1 | (maxangle - minangle) > 5) {
    warning("旋转角度不符合标准，请查询原因", immediate. = T)
  }

  parameter <- list(
    transspot = transspot,
    metaspot = metaspot,
    ratio = meanratio,
    angle = meanangle,
    trans = trans,
    xresolution = xresolution,
    yresolution = yresolution,
    rsdratio = rsdratio,
    rsdangle = rsdangle,
    rawdata = rawdata
  )

  return(parameter)
}

#' @export
getsite <- function(metasite, parameter) {
  metasite2 <- metasite

  # minsite <- which.min((parameter$metaspot[,1]-metasite2[1])^2+(parameter$metaspot[,2]-metasite2[2])^2)
  # minmetasite <- parameter$metaspot[minsite,]
  # minmtranssite <- parameter$transspot[minsite,]

  minmetasite <- parameter$metaspot[1, ]
  minmtranssite <- parameter$transspot[1, ]

  xdistance <- (metasite2[1] - minmetasite[, 1]) * parameter$xresolution / parameter$ratio

  if (parameter$trans) {
    ydistance <- -(metasite2[2] - minmetasite[, 2]) * parameter$yresolution / parameter$ratio
  } else {
    ydistance <- (metasite2[2] - minmetasite[, 2]) * parameter$yresolution / parameter$ratio
  }

  if (ydistance == 0 & xdistance == 0) {
  } else {
    metaangle <- atan2(ydistance, xdistance) / pi * 180
    transangle <- metaangle - parameter$angle
    if (transangle < 0) {
      transangle <- 360 + transangle
    }

    if (transangle < 90 | transangle > 270) {
      xdistance <- sqrt((ydistance^2 + xdistance^2) / (tan(transangle * pi / 180)^2 + 1))
    } else {
      xdistance <- -sqrt((ydistance^2 + xdistance^2) / (tan(transangle * pi / 180)^2 + 1))
    }
    ydistance <- tan(transangle * pi / 180) * xdistance
  }

  xdistance <- (minmtranssite[, 1] / 2 + 0.5) * 100 + xdistance
  ydistance <- minmtranssite[, 2] * sqrt(7500) + ydistance

  transxspot <- (xdistance / 100 - 0.5) * 2
  transyspot <- ydistance / sqrt(7500)

  metaxspot <- xdistance / parameter$xresolution
  metayspot <- ydistance / parameter$yresolution

  return(c(xdistance, ydistance, transxspot, transyspot, metaxspot, metayspot))
}


#' @export
gettranssite <- function(metasite, parameter) {
  metasite2 <- metasite
  metasite2[, c("xdistance", "ydistance", "transxspot", "transyspot", "metaxspot", "metayspot")] <- t(apply(metasite2, 1, getsite, parameter = parameter))
  return(metasite2)
}


#' @export
getalltranssite <- function() {
  numy <- NULL
  for (i in 1:39) {
    numy <- c(numy, rep(c(0, 1), 64) + (i - 1) * 2)
  }
  numx <- rep(0:127, 39)

  transsite <- data.frame(x = numx, y = numy, xdistance = (numx / 2 + 0.5) * 100, ydistance = numy * sqrt(7500))
  return(transsite)
}

#' @export
metatotranssite <- function(transspot,
                            metaspot,
                            xresolution,
                            yresolution,
                            metasite) {
  parameter <- getparameter(
    transspot = transspot,
    metaspot = metaspot,
    xresolution = xresolution,
    yresolution = yresolution
  )
  metasite2 <- gettranssite(
    metasite = metasite,
    parameter = parameter
  )

  data <- list(
    metasite = metasite2,
    parameter = parameter
  )
  return(data)
}

#' imzmltotransdata
#'
#' 将空代与空转数据空间信息配准
#'
#' @param samplename 样本名
#' @param mode 正负离子模式
#' @param xresolution x轴分辨率
#' @param yresolution y轴分辨率
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param map 逻辑，是否保存图片
#' @param coordfile 坐标信息文件
#' @param mzmlpath mzml数据路径
#' @param saveallmzmlpath 转格式后带背景的mzml保存路径
#' @param savemzmlpath 转格式后扣除背景的mzml保存路径
#' @param saveallfinalpath 形成空转分布带背景的mzml保存路径
#' @param savefinalpath 形成空转分布扣除背景的mzml保存路径
#' @param savemappath 转格式后图片保存路径
#' @param savefinalmappath 形成空转分布图片保存路径
#' @param savedatapath 数据矩阵保存路径
#' @param vague mzml数据是否虚化处理
#'
#' @export
imzmltotransdata <- function(samplename,
                             mode = "neg",
                             xresolution = 100,
                             yresolution = 100,
                             mass.range = NULL,
                             resolution = 5,
                             map = T,
                             coordfile = "./sample/union/crood.xlsx",
                             units = "ppm",
                             mzmlpath = "./sample/vague/",
                             saveallmzmlpath = "./sample/union/vague/all/",
                             savemzmlpath = "./sample/union/vague/",
                             saveallfinalpath = "./sample/union/final/all/",
                             savefinalpath = "./sample/union/final/",
                             savemappath = "./sample/union/map/vague/",
                             savefinalmappath = "./sample/union/map/Intensity/",
                             savedatapath = "./sample/union/data/",
                             vague = T) {
  library(Cardinal)
  library(meta)
  
  filename <- paste0(mzmlpath, samplename, "-", mode, ".imzML")
  mse <- lmbio::readdata(filename, attach.only = T,
                         mass.range = mass.range, resolution = resolution, units = units,vague = vague)

  rawdata <- lmbio::readdata(filename = coordfile, sheet = paste0(samplename, "-", mode))
  transspot <- data.frame(x = rawdata$transx, y = rawdata$transy)
  metaspot <- data.frame(x = rawdata$metax, y = rawdata$metay)
  metasite <- data.frame(x = coord(mse)$x, y = coord(mse)$y)

  data <- metatotranssite(transspot = transspot,
                          metaspot = metaspot,
                          xresolution = xresolution,
                          yresolution = yresolution,
                          metasite = metasite)

  coord(mse)$x <- data$metasite$metaxspot
  coord(mse)$y <- data$metasite$metayspot

  if (!dir.exists(saveallmzmlpath)) {
    dir.create(path = saveallmzmlpath, recursive = T)
  }
  filename <- paste0(saveallmzmlpath, samplename, "-", mode, ".imzML")
  filename_ibd <- paste0(saveallmzmlpath, samplename, "-", mode, ".ibd")
  # 判断文件是否存在，存在进行删除
  if (file.exists(filename) | file.exists(filename_ibd)) {
    unlink(filename)
    unlink(filename_ibd)
  }

  coordxy <- data.frame(x = coord(mse)$x,
                        y = coord(mse)$y)

  writeMSIData(mse,
               file = filename,
               outformat = "imzML",
               mz.type = "64-bit float")

  saveRDS(object = coordxy,
          file = paste0(saveallmzmlpath, samplename, "-", mode, ".rds"))


  mse <- mse %>% subsetPixels(pixelApply(mse, sum) != 0)

  filename <- paste0(savemzmlpath, samplename, "-", mode, ".imzML")
  filename_ibd <- paste0(savemzmlpath, samplename, "-", mode, ".ibd")
  # 判断文件是否存在，存在进行删除
  if (file.exists(filename) | file.exists(filename_ibd)) {
    unlink(filename)
    unlink(filename_ibd)
  }

  coordxy <- data.frame(x = coord(mse)$x,
                        y = coord(mse)$y)

  writeMSIData(mse,
               file = filename,
               outformat = "imzML",
               mz.type = "64-bit float")

  saveRDS(object = coordxy,
          file = paste0(savemzmlpath,samplename, "-", mode, ".rds"))

  if (map) {
    Imzmlimage(filename = filename,
               savepath = savemappath,
               savename = paste0("Intensity-", samplename, "-", mode),
               type = c("jpg", "pdf"),
               vague = T,
               vaguerds = paste0(savemzmlpath, samplename, "-", mode, ".rds"))
    Imzmlimage2(filename = filename,
                savepath = paste0(savemappath, samplename, "/", mode),
                savename = paste0("Intensity-", samplename, "-", mode),
                type = c("jpg", "pdf"),
                vague = T,
                vaguerds = paste0(savemzmlpath, samplename, "-", mode, ".rds"))
  }

  transsite <- getalltranssite()

  i <- 1
  disnum <- sqrt((coord(mse)$x * xresolution - transsite[i, "xdistance"])^2 + (coord(mse)$y * yresolution - transsite[i, "ydistance"])^2) < 27.5
  if (any(disnum)) {
    mse4 <- subsetPixels(mse, disnum)
    mse2 <- mse[, 1]
    spectra(mse2)[, ] <- featureApply(mse4, mean)
    coord(mse2)$x <- transsite[i, "xdistance"] / 100
    coord(mse2)$y <- transsite[i, "ydistance"] / 100
    mse2$transx <- transsite[i, "x"]
    mse2$transy <- transsite[i, "y"]
    mse3 <- mse2
    transsite[i, "sample"] <- T
  } else {
    mse2 <- mse[, 1]
    spectra(mse2)[, ] <- 0
    coord(mse2)$x <- transsite[i, "xdistance"] / 100
    coord(mse2)$y <- transsite[i, "ydistance"] / 100
    mse2$transx <- transsite[i, "x"]
    mse2$transy <- transsite[i, "y"]
    mse3 <- mse2
    transsite[i, "sample"] <- F
  }

  for (i in 2:nrow(transsite)) {
    disnum <- sqrt((coord(mse)$x * xresolution - transsite[i, "xdistance"])^2 + (coord(mse)$y * yresolution - transsite[i, "ydistance"])^2) < 27.5
    if (any(disnum)) {
      mse4 <- subsetPixels(mse, disnum)
      mse2 <- mse[, 1]
      spectra(mse2)[, ] <- featureApply(mse4, mean)
      coord(mse2)$x <- transsite[i, "xdistance"] / 100
      coord(mse2)$y <- transsite[i, "ydistance"] / 100
      mse2$transx <- transsite[i, "x"]
      mse2$transy <- transsite[i, "y"]
      mse3 <- BiocGenerics::cbind(mse3, mse2)
      transsite[i, "sample"] <- T
    } else {
      mse2 <- mse[, 1]
      spectra(mse2)[, ] <- 0
      coord(mse2)$x <- transsite[i, "xdistance"] / 100
      coord(mse2)$y <- transsite[i, "ydistance"] / 100
      mse2$transx <- transsite[i, "x"]
      mse2$transy <- transsite[i, "y"]
      mse3 <- BiocGenerics::cbind(mse3, mse2)
      transsite[i, "sample"] <- F
    }
  }

  if (!dir.exists(saveallfinalpath)) {
    dir.create(path = saveallfinalpath, recursive = T)
  }
  filename <- paste0(saveallfinalpath, samplename, "-", mode, ".imzML")
  filename_ibd <- paste0(saveallfinalpath, samplename, "-", mode, ".ibd")
  # 判断文件是否存在，存在进行删除
  if (file.exists(filename) | file.exists(filename_ibd)) {
    unlink(filename)
    unlink(filename_ibd)
  }

  coordxy <- data.frame(x = coord(mse3)$x,
                        y = coord(mse3)$y)

  writeMSIData(mse3,
               file = filename,
               outformat = "imzML",
               mz.type = "64-bit float")

  saveRDS(object = coordxy,
          file = paste0(saveallfinalpath, samplename, "-", mode, ".rds"))

  mse3 <- mse3 %>% subsetPixels(pixelApply(mse3, sum) != 0)

  filename <- paste0(savefinalpath, samplename, "-", mode, ".imzML")
  filename_ibd <- paste0(savefinalpath, samplename, "-", mode, ".ibd")
  # 判断文件是否存在，存在进行删除
  if (file.exists(filename) | file.exists(filename_ibd)) {
    unlink(filename)
    unlink(filename_ibd)
  }

  coordxy <- data.frame(x = coord(mse3)$x,
                        y = coord(mse3)$y)

  writeMSIData(mse3,
               file = filename,
               outformat = "imzML",
               mz.type = "64-bit float")

  saveRDS(object = coordxy,
          file = paste0(savefinalpath, samplename, "-", mode, ".rds"))

  if (!dir.exists(savedatapath)) {
    dir.create(path = savedatapath, recursive = T)
  }

  mz <- data.frame(mz = format(mz(mse3), nsmall = 5))
  transdata <- as.data.frame(iData(mse3, "intensity"))

  names(transdata) <- paste(samplename,
                            mode,
                            mse3$transx, mse3$transy,
                            sep = "-")

  transdata <- cbind(mz, transdata)

  write.table(x = transdata,
              file = paste0(
                savedatapath, samplename,
                "-", mode,
                ".txt"),
              sep = "\t", row.names = F)

  if (map) {
    transimage(filename = filename,
               savepath = savefinalmappath,
               savename = paste0("Intensity-", samplename, "-transmode-", mode),
               type = c("jpg", "pdf"),
               vaguerds = paste0(savefinalpath, samplename, "-", mode, ".rds"))

    transimage2(filename = filename,
                savepath = paste0(savefinalmappath, samplename, "/", mode),
                savename = paste0("Intensity-", samplename, "-transmode-", mode),
                type = c("jpg", "pdf"),
                vaguerds = paste0(savefinalpath, samplename, "-", mode, ".rds"))
  }
}





