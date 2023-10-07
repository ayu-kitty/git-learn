#!/opt/conda/bin/Rscript

#' 通过cluster获取位置
#'
#' @param samplename 样本名称
#' @param mode 正负离子模式
#' @param clusteraera 获取cluster位置
#' @param delclusteraera 删除cluster位置
#' @param clusterpath cluster保存位置
#'
#' @return 位置信息
#'
#' @export
RejectCluster <- function(samplename,
                          mode,
                          clusteraera = NULL,
                          delclusteraera = NULL,
                          clusterpath = "./sample/compress/") {
  ssc <- readRDS(file = paste0(clusterpath,
                               samplename, "-", mode, ".rds"))
  
  ssc <- ssc$sscc
  
  if (!is.null(clusteraera)) {
    samplearea <- ssc$class[[1]] %in% clusteraera
  } else if (!is.null(delclusteraera)) {
    samplearea <- !(ssc$class[[1]] %in% delclusteraera)
  } else {
    samplearea <- ssc$class[[1]] %in% unique(ssc$class[[1]])
  }
  
  return(samplearea)
}


#' 通过强度获取位置
#'
#' @param samplename 样本名称
#' @param mz mz
#' @param limitintensity 限制强度 
#' @param minintensity 最小强度
#' @param imzmlpath 数据路径
#' @param mode 正负离子模式
#'
#' @return 位置信息
#'
#' @export
Rejectintensity <- function(samplename,
                            mode,
                            mz,
                            limitintensity = NULL,
                            minintensity = ifelse(mode == "neg",1000,20000),
                            imzmlpath = "./sample/imzml/") {
  mse <- readimzml(filename = paste0(imzmlpath,"/",samplename,"-",mode,".imzML"))
  mse <- subsetFeatures(x = mse,mz = mz)
  tic <- pixelApply(mse, mean)
  tic2 <- tic
  tic <- smooth.image.gaussian(tic,window=4)
  intensitydata <- tic
  intensitydata <- intensitydata[intensitydata > minintensity]
  
  q3 <- quantile(intensitydata,probs = 0.75)
  q1 <- quantile(intensitydata,probs = 0.25)
  iqr <- q3-q1
  # minlimit <- q1-1.5*iqr
  minlimit <- 0.6*q1
  
  if(!is.null(limitintensity)){
    minlimit <- limitintensity
  }
  
  if(minlimit < minintensity){
    minlimit <- minintensity
  }
  
  print(paste0("使用强度阈值为:",minlimit))
  
  tic[tic < minlimit ] <- 0
  tic[tic >= minlimit ] <- 1
  tic[tic2 < 0.2*minlimit] <- 0
  tic[tic2 >= minlimit] <- 1
  
  img <- image(mse,formula = tic ~ x * y)
  print(img)
  
  samplearea <- as.logical(tic)
  return(samplearea)
}

#' 保存样本位置信息
#'
#' @param samplename 样本名称
#' @param mode 正负离子模式
#' @param clusterpath cluster保存位置
#' @param slidename 玻片名称
#' @param area 区域名称
#' @param savepath 保存位置
#' @param clustertype 图片格式
#' @param ... 见[clusterlimage()]
#'
#' @return 位置信息
#'
#' @export
savearea <- function(slidename,
                     samplename,
                     mode,
                     area,
                     clusterpath = "./sample/compress/",
                     savepath = "./sample/select/",
                     clustertype = "jpg",
                     mirrorx = F,
                     mirrory = F,
                     trans = F,
                     saverotatepath = "./sample/select/rotate/",
                     asp = 1,
                     ...) {
  # 保存area数据
  saverds(data = area,
          filename = paste0(savepath, "/",mode, "/", slidename, "/",samplename, ".rds"))
  
  # 保存旋转数据
  rotatedata <- list(mirrorx = mirrorx,
                     mirrory = mirrory,
                     trans = trans)
  
  saverds(data = rotatedata,
          filename = paste0(saverotatepath, "/",mode, "/", slidename, "/",samplename, ".rds"))
  
  imagearea <<- area
  
  imzmlclusterimage(filename = paste0(clusterpath, slidename, "-", mode, ".rds"),
                    savepath = savepath,
                    mapname = paste0(samplename, "-", mode, "(", slidename, ")"),
                    imagetype = clustertype,
                    subset = imagearea,
                    mirrorx = mirrorx,
                    mirrory = mirrory,
                    trans = trans,
                    asp = asp)
}
