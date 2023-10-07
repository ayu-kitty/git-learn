#!/opt/conda/bin/Rscript

#' 单样本数据平滑处理
#'
#' @param mse 空代数据
#' @param iteration 最大循环次数
#' @param maxPixels 最大点数量
#' @param ...
#'
#' @export
imzmlvague <- function(mse,
                       iteration = 2,
                       maxPixels = 1000000,
                       smooth.image = "gaussian",
                       ...) {
  suppressMessages(library("Cardinal"))
  
  # spectra(mse) <- as.matrix(spectra(mse))
  if(smooth.image != "none"){
    spectra(mse) <- spatialsmooth(mse,window = 3,smooth.image = smooth.image)
  }
  
  mse <- imzmlfill(mse)
  
  for (i in 1:iteration) {
    print(paste0("第",i,"轮虚化"))
    coordxy <- data.frame(
      x = coord(mse)$x,
      y = coord(mse)$y
    )
    intensitydata <- as.matrix(spectra(mse))
    intensitydata <- intensitydata[, order(coordxy$y, coordxy$x)]
    coordxy <- coordxy[order(coordxy$y, coordxy$x), ]
    intensitydata1 <- intensitydata[, coordxy$x != min(coordxy$x)]
    intensitydata2 <- intensitydata[, coordxy$x != max(coordxy$x)]
    intensitydatamean <- (intensitydata1 + intensitydata2) / 2
    
    # for (j in 1:ncol(intensitydatamean)) {
    #   if (sum(intensitydata1[, j]) == 0 & sum(intensitydata2[, j]) != 0) {
    #     intensitydatamean[, j] <- intensitydata1[, j]
    #   } else if (sum(intensitydata1[, j]) != 0 & sum(intensitydata2[, j]) == 0) {
    #     intensitydatamean[, j] <- intensitydata2[, j]
    #   }
    # }
    
    coordxy1 <- coordxy[coordxy$x != min(coordxy$x), ]
    coordxy2 <- coordxy[coordxy$x != max(coordxy$x), ]
    coordxymean <- (coordxy1 + coordxy2) / 2
    mse2 <- mse[, 1:dim(coordxymean)[1]]
    spectra(mse2) <- intensitydatamean
    coord(mse2)$x <- coordxymean$x
    coord(mse2)$y <- coordxymean$y
    
    mse <- BiocGenerics::cbind(mse, mse2)
    
    coordxy <- data.frame(
      x = coord(mse)$x,
      y = coord(mse)$y
    )
    intensitydata <- as.matrix(spectra(mse))
    intensitydata <- intensitydata[, order(coordxy$x, coordxy$y)]
    coordxy <- coordxy[order(coordxy$x, coordxy$y), ]
    intensitydata1 <- intensitydata[, coordxy$y != min(coordxy$y)]
    intensitydata2 <- intensitydata[, coordxy$y != max(coordxy$y)]
    intensitydatamean <- (intensitydata1 + intensitydata2) / 2
    
    # for (j in 1:ncol(intensitydatamean)) {
    #   if (sum(intensitydata1[, j]) == 0 & sum(intensitydata2[, j]) != 0) {
    #     intensitydatamean[, j] <- intensitydata1[, j]
    #   } else if (sum(intensitydata1[, j]) != 0 & sum(intensitydata2[, j]) == 0) {
    #     intensitydatamean[, j] <- intensitydata2[, j]
    #   }
    # }
    
    coordxy1 <- coordxy[coordxy$y != min(coordxy$y), ]
    coordxy2 <- coordxy[coordxy$y != max(coordxy$y), ]
    coordxymean <- (coordxy1 + coordxy2) / 2
    mse2 <- mse[, 1:dim(coordxymean)[1]]
    spectra(mse2) <- intensitydatamean
    coord(mse2)$x <- coordxymean$x
    coord(mse2)$y <- coordxymean$y
    
    mse <- BiocGenerics::cbind(mse, mse2)
    if (dim(mse)[2] > maxPixels) {
      break
    }
  }
  
  return(mse)
}