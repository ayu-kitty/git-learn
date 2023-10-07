#!/opt/conda/bin/Rscript
#' getfigdataforresnet
#'
#' 自动加手动调整HE染色图
#'
#' @param imzmlfile imzml数据名称
#' @param samplename 样本名称
#' @param mode 正负离子模式
#' @param figrds 图像名称
#' @param figname 图片名称
#'
#' @export
get_fig_data_for_resnet <- function(imzmlfile,
                                mode = ifelse(grepl(pattern = "-neg\\.",basename(imzmlfile)),"neg","pos"),
                                samplename = gsub(pattern = paste0("-",mode,"\\..*$"),replacement = "",x = basename(imzmlfile)),
                                figrds = NULL,
                                figname = NULL){
  library(EBImage)
  library(lmbio)

  mse <- readimzml(filename = imzmlfile)

  spectradata <- as.matrix(iData(mse, "intensity"))
  spectradata <- as.data.frame(spectradata)
  colnames(spectradata) <- paste(coord(mse)$x, coord(mse)$y,
                                 sep = "-")
  row.names(spectradata) <- format(mz(mse), nsmall = 5, trim = T)
  infodata <- data.frame(Name = colnames(spectradata),
                     Samplename = samplename,
                     mode = mode,
                     x = coord(mse)$x,
                     y = coord(mse)$y,
                     stringsAsFactors = F)

  spectradata <- spectradata[!apply(spectradata,MARGIN = 1,FUN = function(x){all(x==0)}),]
  # 1. 数据平滑
  spectradata <- t(apply(spectradata,MARGIN = 1,smooth.image.gaussian))

  # 2. 数据标准化
  spectradata <- t(scale(t(spectradata),center = T,scale = T))
  spectradata[is.na(spectradata)] <- 0

  # 3. 数据归一化
  spectradata2 <- apply(X = spectradata,MARGIN = 1,scalezerotoone)
  spectradata <- t(spectradata2)

  spectradata <- as.data.frame(spectradata)
  colnames(spectradata) <- paste(coord(mse)$x, coord(mse)$y,
                                 sep = "-")

  resnetdata <- list(spectradata = spectradata,
                     infodata = infodata,
                     anadata = list(r = spectradata,
                                    g = spectradata,
                                    b = spectradata))

  # 4. 获取figrds
  if(!is.null(figrds)){
    figresizedata_r_1 <- infodata[,c("x","y")]
    figresizedata_g_1 <- infodata[,c("x","y")]
    figresizedata_b_1 <- infodata[,c("x","y")]

    xlength <- max(coord(mse)$x)-min(coord(mse)$x)+1
    ylength <- max(coord(mse)$y)-min(coord(mse)$y)+1

    for ( i in 1:length(figrds)) {
      figdata <- readRDS(figrds[i])
      figresizedata <- EBImage::resize(x = figdata$figcut,w = xlength,h = ylength)

      figresizedata <- imageData(figresizedata)

      figresizedata_r <- t(figresizedata[,,1])
      colnames(figresizedata_r) <- min(coord(mse)$x):max(coord(mse)$x)
      rownames(figresizedata_r) <- min(coord(mse)$y):max(coord(mse)$y)
      figresizedata_r <- reshape2::melt(data = figresizedata_r,varnames=c("y","x"),value.name = paste0(figname[i]))
      figresizedata_r <- figresizedata_r[paste0(figresizedata_r$x,"-",figresizedata_r$y) %in% paste0(infodata$x,"-",infodata$y),]
      figresizedata_r_1 <- merge(figresizedata_r_1,figresizedata_r,by = c("x","y"))

      figresizedata_g <- t(figresizedata[,,2])
      colnames(figresizedata_g) <- min(coord(mse)$x):max(coord(mse)$x)
      rownames(figresizedata_g) <- min(coord(mse)$y):max(coord(mse)$y)
      figresizedata_g <- reshape2::melt(data = figresizedata_g,varnames=c("y","x"),value.name = paste0(figname[i]))
      figresizedata_g <- figresizedata_g[paste0(figresizedata_g$x,"-",figresizedata_g$y) %in% paste0(infodata$x,"-",infodata$y),]
      figresizedata_g_1 <- merge(figresizedata_g_1,figresizedata_g,by = c("x","y"))

      figresizedata_b <- t(figresizedata[,,3])
      colnames(figresizedata_b) <- min(coord(mse)$x):max(coord(mse)$x)
      rownames(figresizedata_b) <- min(coord(mse)$y):max(coord(mse)$y)
      figresizedata_b <- reshape2::melt(data = figresizedata_b,varnames=c("y","x"),value.name = paste0(figname[i]))
      figresizedata_b <- figresizedata_b[paste0(figresizedata_b$x,"-",figresizedata_b$y) %in% paste0(infodata$x,"-",infodata$y),]
      figresizedata_b_1 <- merge(figresizedata_b_1,figresizedata_b,by = c("x","y"))
    }

    figresizedata_r_1_1 <- t(figresizedata_r_1[,-1:-2,drop =F])
    colnames(figresizedata_r_1_1) <- paste0(figresizedata_r_1$x,"-",figresizedata_r_1$y)
    figresizedata_r_1_1 <- figresizedata_r_1_1[,colnames(spectradata),drop=F]
    spectradata_r <- rbind(spectradata,figresizedata_r_1_1)

    figresizedata_g_1_1 <- t(figresizedata_g_1[,-1:-2,drop =F])
    colnames(figresizedata_g_1_1) <- paste0(figresizedata_g_1$x,"-",figresizedata_g_1$y)
    figresizedata_g_1_1 <- figresizedata_g_1_1[,colnames(spectradata),drop=F]
    spectradata_g <- rbind(spectradata,figresizedata_g_1_1)

    figresizedata_b_1_1 <- t(figresizedata_b_1[,-1:-2,drop =F])
    colnames(figresizedata_b_1_1) <- paste0(figresizedata_b_1$x,"-",figresizedata_b_1$y)
    figresizedata_b_1_1 <- figresizedata_b_1_1[,colnames(spectradata),drop=F]
    spectradata_b <- rbind(spectradata,figresizedata_b_1_1)

    resnetdata$figdata$r <- figresizedata_r_1_1
    resnetdata$figdata$g <- figresizedata_g_1_1
    resnetdata$figdata$b <- figresizedata_b_1_1

    resnetdata$anadata$r <- spectradata_r
    resnetdata$anadata$g <- spectradata_g
    resnetdata$anadata$b <- spectradata_b
  }

  return(resnetdata)
}

#' @export
scalezerotoone <- function(x){
  (x - min(x))/(max(x)-min(x))
}

#' @export
smooth.image.gaussian <- function(x, window=3, ...) {
  if ( all(is.na(x)) ) return(x)
  r <- floor(window / 2)
  sd <- window / 4
  x.new <- .Call("C_gaussianFilter", x, r, sd, PACKAGE="Cardinal")
  x.new <- max(x, na.rm=TRUE) * x.new / max(x.new, na.rm=TRUE)
  x.new
}

