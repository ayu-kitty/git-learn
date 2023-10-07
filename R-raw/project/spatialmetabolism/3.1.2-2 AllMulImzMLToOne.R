#!/opt/conda/bin/Rscript

#' 将目录下HF RAW 数据按文件目录进行多imzml合并为一
#'
#' @param imzmlpath 数据目录
#' @param moderange 正负离子模式
#' @param ... 见[MulImzMLToOne()]
#'
#' @export
AllMulImzMLToOne <- function(imzmlpath = "./raw/",
                             moderange = c("neg", "pos"),
                             ratio = 2,
                             ...) {
  for (mode in moderange) {
    
    # 获取.imzML文件绝对路径
    filename <- list.files(path = paste0(imzmlpath, mode),
                           full.names = T)
    
    if (length(filename) == 0) {
      print(paste0("在",imzmlpath,"目录下",mode,"未找到目录"))
      next
    }
    
    # 获取玻片名
    slidename <- list.files(path = paste0(imzmlpath,"/",mode),
                            full.names = F)
    
    for (i in 1:length(filename)) {
      ImzMLToOne(imzMLdir = filename[i],
                 slidename = slidename[i],
                 mode = mode,
                 ratio = ratio,
                 ...)
    }
  }
}

#' 多imzml合并为一
#'
#' @param imzMLdir 文件路径
#' @param slidename 玻片名
#' @param mode 正负离子模式
#' @param datatrans 是否转置
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param savepath 保存路径
#' @param imagetype 图片格式
#' @param ... 见[Imzmlimage()]
#'
#' @export
ImzMLToOne <- function(imzMLdir,
                       slidename,
                       mode,
                       datatrans = F,
                       mass.range = NULL,
                       resolution = 5,
                       units = "ppm",
                       savepath = "./sample/imzml/",
                       imagetype = "jpg",
                       ratio = 2,
                       asp = 1,
                       ...) {
  suppressMessages(library("Cardinal"))
  print(paste0(slidename,"运行中"))
  
  # 读取数据量list
  imzMLfile <- list.files(path = imzMLdir, pattern = ".imzML$", full.names = T)
  imzMLfile2 <- list.files(path = imzMLdir, pattern = ".imzML$", full.names = F)
  mzMLfile <- gsub(pattern = "imzML$", replacement = "mzML", x = imzMLfile)
  
  if (!all(file.exists(mzMLfile))) {
    stop("mzML数据不全")
  }
  
  # 判断目录下是否有imzML文件
  if (length(imzMLfile) == 0) {
    return(paste0(slidename, "不存在imzML数据"))
  } else {
    
    # imzML文件顺序是否装置
    if (datatrans) {
      imzMLfile <- imzMLfile[length(imzMLfile):1]
      imzMLfile2 <- imzMLfile2[length(imzMLfile2):1]
      mzMLfile <- mzMLfile[length(mzMLfile):1]
    }
    
    # 第一个imzML读取
    mse <- readMSIData(file = imzMLfile[1], 
                       attach.only = F, 
                       resolution = resolution, 
                       units = units,
                       mass.range = mass.range)
    
    # 提取resolution值，比如如果默认为5ppm，则此处提取的resolution值就是5；除非如果人为没定义resolution值，即NA，程序则会自动从数据中猜测一个分辨率。
    resolution <- resolution(mse)
    
    centroided(mse) <- F
    
    # 读取ms1光谱
    mzml <- MSnbase::readMSData(files = mzMLfile[1], msLevel. = 1)
    
    # 用保留时间换算像素点位次round
    rttime <- data.frame(rttime = rtime(mzml) * ratio)
    rttime[, "rttime"] <- rttime$rttime - rttime$rttime[1] + 1
    rttime[, "round"] <- round(rttime$rttime)
    
    # 把位次转换为x坐标
    coord(mse)$x <- rttime[, "round"]
    
    # 计算长度
    range <- min(coord(mse)$x):max(coord(mse)$x)
    
    # 缺位和多余像素点提出
    missingrange <- range[!(range %in% coord(mse)$x)]
    
    # 空缺像素点处理,多余像素点删除,空缺像素点补充(*)
    if (length(missingrange) != 0) {
      addrange <- NULL
      for (k in seq_len(length(missingrange))) {
        addrange <- c(addrange, which.min(abs(coord(mse)$x - missingrange[k])))
      }
      addmse <- mse[, addrange]
      coord(addmse)$x <- missingrange
      mse <- BiocGenerics::cbind(mse, addmse)
    }
    
    # 去重
    mse <- mse[, !duplicated(coord(mse)$x)]
    
    
    # imzML读取读取合并，并根据读取顺序对y赋值,其他同上
    if (length(imzMLfile) > 1) {
      for (j in 2:length(imzMLfile)) {
        mse2 <- readMSIData(imzMLfile[j], 
                            attach.only = F, 
                            resolution = resolution, 
                            units = units,
                            mass.range = mass.range)
        coord(mse2)$y <- j
        centroided(mse2) <- F
        
        mzml <- MSnbase::readMSData(files = mzMLfile[j], msLevel. = 1)
        rttime <- data.frame(rttime = rtime(mzml) * ratio)
        rttime[, "rttime"] <- rttime$rttime - rttime$rttime[1] + 1
        rttime[, "round"] <- round(rttime$rttime)
        coord(mse2)$x <- rttime[, "round"]
        range <- min(coord(mse2)$x):max(coord(mse2)$x)
        missingrange <- range[!(range %in% coord(mse2)$x)]
        if (length(missingrange) != 0) {
          addrange <- NULL
          for (k in seq_len(length(missingrange))) {
            addrange <- c(addrange, which.min(abs(coord(mse2)$x - missingrange[k])))
          }
          addmse <- mse2[, addrange]
          coord(addmse)$x <- missingrange
          mse2 <- BiocGenerics::cbind(mse2, addmse)
        }
        mse2 <- mse2[, !duplicated(coord(mse2)$x)]
        mse <- BiocGenerics::cbind(mse, mse2)
      }
    }
    
    # 数据整理合并
    mse <- mse[, order(coord(mse)$y, coord(mse)$x)]
    
    coordxy <- coord(mse)
    
    coordxy <- data.frame(
      x = coordxy$x,
      y = coordxy$y
    )
    coordxysum <- NULL
    
    for (f in unique(coordxy$y)) {
      maxx <- max(coordxy[coordxy$y == f, "x"])
      coordxysum2 <- data.frame(x = maxx, y = f)
      coordxysum <- rbind(coordxysum, coordxysum2)
    }
    
    for (j in seq_len(nrow(coordxysum))) {
      while (coordxysum$x[j] != max(coordxysum$x)) {
        print(coordxysum$y[j])
        mse2 <- mse[, (coordxy$x == coordxysum$x[j]) & (coordxy$y == coordxysum$y[j])]
        coord(mse2)$x <- coordxysum$x[j] + 1
        mse <- BiocGenerics::cbind(mse, mse2)
        mse <- mse[, order(coord(mse)$y, coord(mse)$x)]
        coordxy <- coord(mse)
        coordxy <- data.frame(
          x = coordxy$x,
          y = coordxy$y
        )
        coordxysum <- NULL
        for (f in unique(coordxy$y)) {
          maxx <- max(coordxy[coordxy$y == f, "x"])
          coordxysum2 <- data.frame(x = maxx, y = f)
          coordxysum <- rbind(coordxysum, coordxysum2)
        }
      }
    }
    
    # 运行集赋予名称
    run(mse) <- slidename
    
    filename <- paste0(savepath,"/",slidename,"-",mode,".imzML")
    
    saveimzml(data = mse,
              filename = filename)
    
    imzmlimage(filename = mse,
               savepath = savepath,
               mapname = paste0(slidename, "-", mode),
               imagetype = imagetype,
               asp = asp)
    
    getimzmlbasicinfo(mse = mse,
                      mode = mode,
                      filename = paste0(savepath,"/",slidename, "-", mode,".xlsx"))
    getimzmlbasicinfo(mse = mse[,1:10],
                      mode = mode,
                      filename = paste0(savepath,"/",slidename, "-", mode,"-bg.xlsx"))
    
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-i","--imzmlpath",default = "./raw/", help = "imzml原始文件路径,默认./raw/")
  parser$add_argument("-s","--savepath",default = "./sample/imzml/", help = "imzml数据保存路径,默认./sample/imzml/")
  parser$add_argument("-r","--ratio",default = 2, type= "double",help = "仪器每秒扫描数量,默认2",required = T)
  parser$add_argument("-m","--massrange",default = NULL,type= "double",help = "mz范围,默认NULL",
                      nargs = 2,dest = "mass.range")
  parser$add_argument("-re","--resolution",default = 5, type= "double",help = "分辨率,默认5")
  parser$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  parser$add_argument("-a","--asp",default = 1, type= "double",help = "长宽比，默认为1")
  args <- parser$parse_args()
  
  writeinfo()
  createdir(filename = args$savepath,linkdir = T)
  
  mulargs <- do.call(what = AllMulImzMLToOne,args = args)
  
  writeinfo(endtime = T)
  
}
