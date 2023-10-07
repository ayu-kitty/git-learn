#!/opt/conda/bin/Rscript

#' sscc运算
#'
#' @param imzmlpath 数据路径
#' @param area 是否绘制选择区域
#' @param areaname 区域名
#' @param areapath 区域信息文件路径
#' @param r 距离半径
#' @param k 聚类数量
#' @param s 稀疏参数
#' @param iter.max 聚类循环次数
#' @param samplename 样本名
#' @param mode 正负离子模式
#' @param savename 文件保存名字
#' @param savepath 文件保存路径
#' @param filtermz 根据提供的mz运算
#' @param ... 见[spatialShrunkenCentroids()]
#'
#' @export
imzmltosscc <- function(imzmlpath = "./sample/final/",
                        samplename = NULL,
                        savename = samplename[1],
                        mode = "neg",
                        r = 1,
                        k = 8,
                        s = 3,
                        iter.max = 30,
                        savepath = "./sample/cluster/sscc/",
                        area = F,
                        areaname = "ALL",
                        areapath = "./area/data/",
                        filtermz = NULL,
                        smooth.image = "none",
                        filename = paste0(imzmlpath, samplename, "-", mode, ".imzML"),
                        areards = paste0(areapath,samplename,"-",areaname,"-", mode, ".rds"),
                        savedata = T,
                        ...) {
  suppressMessages(library("Cardinal"))
  setCardinalBPPARAM(MulticoreParam(workers = 3))
  
  # filename <- paste0(imzmlpath, samplename, "-", mode, ".imzML")
  # areards <- paste0(areapath,samplename,"-",areaname,"-", mode, ".rds")
  
  #读取样本imzml
  mse <- readdata(filename = filename,
                  area = area,
                  areards = areards,
                  filtermz = filtermz,
                  attach.only = F,
                  addy = T)
  
  #平滑化处理
  if(smooth.image != "none"){
    spectra(mse) <- spatialsmooth(mse,window = 2,smooth.image = smooth.image)
  }
  
  print("数据提取中")
  #提取表达量
  spectradata <- as.matrix(iData(mse, "intensity"))
  spectradata <- as.data.frame(spectradata)
  #命名 样本名-mode-x坐标-y坐标
  colnames(spectradata) <- paste(mse$samplename, mode,
                                 coord(mse)$x, coord(mse)$y,
                                 sep = "-")
  #mz保留5位小数
  row.names(spectradata) <- format(mz(mse), nsmall = 5, trim = T)
  rawspectradata <- spectradata
  # spectradata <- t(scale(t(spectradata),center = T,scale = T))
  # spectradata[is.na(spectradata)] <- 0
  # spectra(mse) <- spectradata
  
  print("sscc聚类分析开始")
  #空间收缩质心  r:附近像素的空间邻域半径  s:稀疏性阈值参数，用于收缩t统计量 iter.max:聚类循环次数
  ssc <- spatialShrunkenCentroids(mse,
                                  r = r,
                                  k = k,
                                  s = s,
                                  iter.max = iter.max,
                                  ...)
  print("sscc聚类分析结束")
  
  #数据格式：聚类数/x坐标/y坐标/样本名/mode/像素点名称
  plotdata <- data.frame(Cluster = ssc@resultData@listData[[1]][["class"]],
                         x = coord(mse)$x,
                         y = coord(mse)$y,
                         Sample = run(mse),
                         mode = mode,
                         Name = colnames(spectradata),
                         stringsAsFactors = F)

  alldata <- list(spectradata = rawspectradata,
                  sscc = ssc,
                  plotdata = plotdata)
  
  setCardinalBPPARAM(SerialParam())
  if(savedata){
    # 数据保存
    saverds(data = alldata,
            filename = paste0(savepath,"/",savename, "-", mode, "-data.rds"))
  }else{
    return(alldata)
  }

}

#' tsne运算
#'
#' @param imzmlpath 数据路径
#' @param area 是否绘制选择区域
#' @param areaname 区域名
#' @param areapath 区域信息文件路径
#' @param samplename 样本名
#' @param mode 正负离子模式
#' @param savename 文件保存名字
#' @param savepath 文件保存路径
#' @param filtermz 根据提供的mz运算
#' @param ... 见[Rtsne()]
#'
#' @export
imzmltotsne <- function(imzmlpath = "./sample/final/",
                        samplename = NULL,
                        savename = samplename[1],
                        mode = "neg",
                        savepath = "./sample/cluster/tsne/",
                        area = F,
                        areaname = "ALL",
                        areapath = "./area/data/",
                        dims = 2,
                        verbose = F,
                        pca = T,
                        check_duplicates = F,
                        filtermz = NULL,
                        smooth.image = "none",
                        ...) {
  suppressMessages(library("Cardinal"))
  suppressMessages(library("Rtsne"))
  
  filename <- paste0(imzmlpath, samplename, "-", mode, ".imzML")
  areards <- paste0(areapath,samplename,"-",areaname,"-", mode, ".rds")
  
  mse <- readimzml(filename = filename,
                   area = area,
                   areards = areards,
                   filtermz = filtermz,
                   attach.only = F,
                   addy = T)
  
  #平滑化处理
  if(smooth.image != "none"){
    spectra(mse) <- spatialsmooth(mse,window = 2,smooth.image = smooth.image)
  }
  
  print("数据提取中")
  #提取表达量
  spectradata <- as.matrix(iData(mse, "intensity"))
  
  spectradata <- as.data.frame(spectradata)
  #命名 样本名-mode-x坐标-y坐标
  colnames(spectradata) <- paste(mse$samplename, mode,
                                 coord(mse)$x, coord(mse)$y,
                                 sep = "-")
  #mz保留5位小数
  row.names(spectradata) <- format(mz(mse), nsmall = 5, trim = T)
  rawspectradata <- spectradata
  caldata <- t(spectradata)
  print(length(caldata))
  print(dim(caldata))
  
  # if(length(caldata) > 500000000){
  #   print(dim(caldata))
  #   print("进行pca降维运算")
  #   pca1 <- prcomp(caldata,center = TRUE,scale. = F,rank. = 50)
  #   caldata <- as.matrix(pca1$x)
  #   pca <- F
  #   print(dim(caldata))
  # }
  
  print("进行pca降维运算")
  # caldata <<- caldata
  # caldata <- caldata[,!(apply(caldata,MARGIN = 2,FUN = sd)==0)]
  # pca1 <- prcomp(caldata,center = TRUE,scale. = TRUE,rank. = 20)
  # caldata <- as.matrix(pca1$x)
  print("tsne运算中")
  #dims:输出维度  initial_dims:初始PCA步骤中保留的维度数  verbose:是否打印进度更新 pca:是否执行初始PCA步骤  check_duplicates:检查是否存在重复项 perplexity:衡量概率模型预测能力的指标，控制降维后数据点之间的距离
  tsne_i <- Rtsne(caldata,
                  dims = dims,
                  initial_dims = 50,
                  verbose = verbose,
                  pca = pca,
                  check_duplicates = check_duplicates,
                  perplexity = ifelse(dim(caldata)[1] < 91,10,30))
  print("tsne运算结束")
  #stringsAsFactors:在读入数据遇到字符串时，不将其转换为factors，仍然保留为字符串格式
  #数据格式：TSNE维度坐标/x坐标/y坐标/样本名/mode/像素点名称
  plotdata <- data.frame(tsne_i$Y,
                         x = coord(mse)$x,
                         y = coord(mse)$y,
                         Sample = run(mse),
                         mode = mode,
                         Name = colnames(spectradata),
                         stringsAsFactors = F)
  names(plotdata)[1:dims] <- paste("Tsne", 1:dims)
  
  alldata <- list(spectradata = rawspectradata,
                  tsne = tsne_i,
                  plotdata = plotdata)
  
  # 数据保存
  print("tsne降维数据保存中")
  saverds(data = alldata,
          filename = paste0(savepath,"/",savename, "-", mode, "-data.rds"))
}

#' umap运算
#'
#' @param imzmlpath 数据路径
#' @param area 是否绘制选择区域
#' @param areaname 区域名
#' @param areapath 区域信息文件路径
#' @param samplename 样本名
#' @param mode 正负离子模式
#' @param smoothSignal_method 平滑模式
#' @param savename 文件保存名字
#' @param savepath 文件保存路径
#' @param filtermz 根据提供的mz运算
#' @param ... 见[umap()]
#'
#' @export
imzmltoumap <- function(imzmlpath = "./sample/final/",
                        samplename = NULL,
                        savename = samplename[1],
                        mode = "neg",
                        smoothSignal_method = NULL,
                        savepath = "./sample/cluster/umap/",
                        area = F,
                        areaname = "ALL",
                        areapath = "./area/data/",
                        n_components = 2,
                        filtermz = NULL,
                        smooth.image = "none",
                        ...) {
  suppressMessages(library("Cardinal"))
  suppressMessages(library("umap"))
  
  filename <- paste0(imzmlpath, samplename, "-", mode, ".imzML")
  areards <- paste0(areapath,samplename,"-",areaname,"-", mode, ".rds")
  
  mse <- readimzml(filename = filename,
                   area = area,
                   areards = areards,
                   filtermz = filtermz,
                   attach.only = F,
                   addy = T)
  #平滑化处理
  if(smooth.image != "none"){
    spectra(mse) <- spatialsmooth(mse,window = 2,smooth.image = smooth.image)
  }
  
  print("数据提取中")
  #提取表达量
  spectradata <- as.matrix(iData(mse, "intensity"))
  spectradata <- as.data.frame(spectradata)
  #命名 样本名-mode-x坐标-y坐标
  colnames(spectradata) <- paste(mse$samplename, mode,
                                 coord(mse)$x, coord(mse)$y,
                                 sep = "-")
  #mz保留5位小数
  row.names(spectradata) <- format(mz(mse), nsmall = 5, trim = T)
  rawspectradata <- spectradata
  caldata <- t(spectradata)
  print(length(caldata))
  
  # if(length(caldata) > 500000000){
  #   print("进行pca降维运算")
  #   pca1 <- prcomp(caldata,center = TRUE,scale. = TRUE,rank. = 50)
  #   caldata <- as.matrix(pca1$x)
  # }
  
  print("进行pca降维运算")
  # caldata <- caldata[,!(apply(caldata,MARGIN = 2,FUN = sd)==0)]
  # pca1 <- prcomp(caldata,center = TRUE,scale. = TRUE,rank. = 20)
  # caldata <- as.matrix(pca1$x) 
  print("umap运算中")
  #n_components:空间维度  init:坐标的初始化类型  pca:将数据减少到此列数 pca_center:计算PCA之前将数据中心化
  umap_i <- uwot::umap(caldata,
                       n_components = n_components,
                       # init = "random",
                       init = "spectral",
                       pca = 50,
                       pca_center = TRUE)
  print("umap运算结束")
  #数据格式：UMAP维度坐标/x坐标/y坐标/样本名/mode/像素点名称
  plotdata <- data.frame(umap_i,
                         x = coord(mse)$x,
                         y = coord(mse)$y,
                         Sample = run(mse),
                         mode = mode,
                         Name = colnames(spectradata),
                         stringsAsFactors = F)
  names(plotdata)[1:n_components] <- paste("Umap", 1:n_components)
  
  alldata <- list(spectradata = rawspectradata,
                  umap = umap_i,
                  plotdata = plotdata)
  
  # 数据保存
  print("umap降维数据保存中")
  saverds(data = alldata,
          filename = paste0(savepath,"/",savename, "-", mode, "-data.rds"))
}


#' 空代单样本聚类与降维分析
#'
#' @param imzmlpath 数据路径
#' @param area 是否绘制选择区域
#' @param areaname 区域名
#' @param areapath 区域信息文件路径
#' @param moderange 正负离子模式
#' @param moudle 运行模块
#' @param ... 见moudle函数
#'
#' @export
MulsingleCluster <- function(imzmlpath = "./sample/final/",
                             area = F,
                             areaname = "ALL",
                             areapath = "./sample/area/data/",
                             moderange = c("neg", "pos"),
                             moudle = imzmltosscc,
                             wantsample = NULL,
                             mulname = NULL,
                             qcanalyst = F,
                             ...) {
  for (mode in moderange) {
    
    # 获取.imzML文件绝对路径
    filename <- list.files(path = imzmlpath,
                           pattern = paste0("-", mode, ".imzML$"),
                           full.names = F,
                           recursive = T)
    
    if(length(filename) == 0){
      print(paste0("在",imzmlpath,"目录下未找到", mode, "模式的imzML文件"))
      next
    }
    
    # 获取样品名
    samplename <- gsub(pattern = paste0("-", mode, ".imzML"),
                       replacement = "",
                       x = filename)
    if(!is.null(wantsample)){
      samplename <- samplename[samplename %in% wantsample]
    }
    
	#如果不做QC分析则剔除QC样本名
    if(!qcanalyst){
      delsample2 <- "qc_data"
      if(any(grepl(pattern = paste0("^",delsample2,".*"),x = samplename))){
        samplename <- samplename[!grepl(pattern = paste0("^",delsample2,".*"),x = samplename)]
      }
    }
    
    imagename <- samplename
    
    areaname2 <- NULL
    samplename2 <- NULL
	#是否进行选区的聚类分析
    if(area){
      for ( i in seq_len(length(samplename))) {
        areards <- paste0(areapath,samplename[i],"-",areaname,"-", mode, ".rds")
        areaname3 <- areaname[file.exists(areards)]
        samplename3 <- rep(samplename[i],length(areaname3))
        areaname2 <- c(areaname2,areaname3)
        samplename2 <- c(samplename2,samplename3)
      }
      samplename <- samplename2
      imagename <- paste0(samplename,"-",areaname2)
    }
    
    # 判断是否找到imzML文件
    if (length(samplename) > 0) {
      
      
      # 针对imzML文件循环剔除背景处理
      for (i in seq_len(length(samplename))) {
        
        ssc <- moudle(samplename = samplename[i],
                      savename = imagename[i],
                      imzmlpath = imzmlpath,
                      mode = mode,
                      area = area,
                      areaname = areaname2[i],
                      areapath =  areapath,
                      ...)
        gc(reset = TRUE)
      }
    } else {
      print(paste0("在",imzmlpath,"目录下未找到", mode, "模式的imzML文件"))
    }
  }
}

#' 空代单样本聚类与降维分析
#'
#' @param imzmlpath 数据路径
#' @param area 是否绘制选择区域
#' @param areaname 区域名
#' @param areapath 区域信息文件路径
#' @param moderange 正负离子模式
#' @param moudle 运行模块
#' @param mulname 多样本名称
#' @param ... 见moudle函数
#'
#' @export
MulallCluster <- function(imzmlpath = "./sample/final/",
                          area = F,
                          areaname = "ALL",
                          areapath = "./sample/area/data/",
                          moderange = c("neg", "pos"),
                          moudle = imzmltosscc,
                          mulname = ifelse(area,paste0("all-",areaname[1]),"all"),
                          k = 15,
                          wantsample = NULL,
                          qcanalyst = F,
                          ...) {
  for (mode in moderange) {
    
    # 获取.imzML文件绝对路径
    filename <- list.files(path = imzmlpath,
                           pattern = paste0("-", mode, ".imzML$"),
                           full.names = F,
                           recursive = T)
    
    if(length(filename) == 0){
      print(paste0("在",imzmlpath,"目录下未找到", mode, "模式的imzML文件"))
      next
    }
    
    # 获取样品名
    samplename <- gsub(pattern = paste0("-", mode, ".imzML"),
                       replacement = "",
                       x = filename)
    
    if(!is.null(wantsample)){
      samplename <- wantsample
    }
    #如果不做QC分析则剔除QC样本名
    delsample2 <- "qc_data"
    if(qcanalyst){
      if(any(grepl(pattern = paste0("^",delsample2,".*"),x = samplename))){
      }else{
        return()
      }
    }else{
      if(any(!grepl(pattern = paste0("^",delsample2,".*"),x = samplename))){
        samplename <- samplename[!grepl(pattern = paste0("^",delsample2,".*"),x = samplename)]
      }else{
        return()
      }
    }
    
    imagename <- samplename
    
    areaname2 <- NULL
    samplename2 <- NULL
	#是否进行选区的聚类分析
    if(area){
      for ( i in seq_len(length(samplename))) {
        areards <- paste0(areapath,samplename[i],"-",areaname,"-", mode, ".rds")
        areaname3 <- areaname[file.exists(areards)]
        samplename3 <- rep(samplename[i],length(areaname3))
        areaname2 <- c(areaname2,areaname3)
        samplename2 <- c(samplename2,samplename3)
      }
      samplename <- samplename2
      imagename <- paste0(samplename,"-",areaname2)
    }
    
    # 判断是否找到imzML文件
    if (length(samplename) > 1) {
      #根据moudle传参做对应分析
      ssc <- moudle(samplename = samplename,
                    savename = mulname,
                    imzmlpath = imzmlpath,
                    mode = mode,
                    area = area,
                    areaname = areaname2,
                    areapath =  areapath,
                    k = k,
                    ...)
      gc(reset = TRUE)
      
    } else {
      print(paste0("在",imzmlpath,"目录下未找到", mode, "模式的imzML文件"))
    }
  }
}

#' 进行sscc、umap、tsne运算
#'
#' @param imzmlpath 数据路径，默认./sample/adjustdata/
#' @param savepath 保存路径
#' @param ssccsavepath sscc保存路径，默认./sample/cluster/sscc/
#' @param tsnesavepath tsne保存路径，默认./sample/cluster/tsne/
#' @param umapsavepath umap保存路径，默认./sample/cluster/umap/
#' @param k 聚类数，默认8
#' @param mulk 多组聚类数，默认15
#' @param singleanalyst 是否单样本分析
#' @param allanalyst 是否多样本分析
#' @param ...
#'
#' @export
MulanaCluster <- function(imzmlpath = "./sample/adjustdata/",
                          savepath = "./sample/cluster/",
                          ssccsavepath = paste0(savepath,"sscc/"),
                          tsnesavepath = paste0(savepath,"tsne/"),
                          umapsavepath = paste0(savepath,"umap/"),
                          k = 8,
                          s = 3,
                          r = 1,
                          mulk = 15,
                          muls = 3,
                          mulr = 1,
                          singleanalyst = T,
                          allanalyst = T,
                          ...){
  if(singleanalyst){
    # 单样本分析
	#sscc分析
    MulsingleCluster(moudle = imzmltosscc,
                     imzmlpath = imzmlpath,
                     k = k,
                     s = s,
                     r = r,
                     savepath = ssccsavepath,
                     ...)
    #tsne分析
	MulsingleCluster(moudle = imzmltotsne,
                     imzmlpath = imzmlpath,
                     savepath = tsnesavepath,
                     smooth.image = "gaussian",
                     ...)
    #umap分析
	MulsingleCluster(moudle = imzmltoumap,
                     imzmlpath = imzmlpath,
                     savepath = umapsavepath,
                     smooth.image = "gaussian",
                     ...)
  }
  
  if(allanalyst){
    # 多样本分析
    MulallCluster(moudle = imzmltosscc,
                  imzmlpath = imzmlpath,
                  k = mulk,
                  s = muls,
                  r = mulr,
                  savepath = ssccsavepath,
                  smooth.image = "gaussian",
                  ...)
    MulallCluster(moudle = imzmltotsne,
                  imzmlpath = imzmlpath,
                  savepath = tsnesavepath,
                  smooth.image = "gaussian",
                  ...)
    MulallCluster(moudle = imzmltoumap,
                  imzmlpath = imzmlpath,
                  savepath = umapsavepath,
                  smooth.image = "gaussian",
                  ...)
    
    # 带QC分析
    MulallCluster(moudle = imzmltosscc,
                  imzmlpath = imzmlpath,
                  k = mulk,
                  s = muls,
                  r = mulr,
                  savepath = ssccsavepath,
                  smooth.image = "gaussian",
                  qcanalyst = T,
                  mulname = "all-qc",
                  ...)
    MulallCluster(moudle = imzmltotsne,
                  imzmlpath = imzmlpath,
                  savepath = tsnesavepath,
                  smooth.image = "gaussian",
                  qcanalyst = T,
                  mulname = "all-qc",
                  ...)
    MulallCluster(moudle = imzmltoumap,
                  imzmlpath = imzmlpath,
                  savepath = umapsavepath,
                  smooth.image = "gaussian",
                  qcanalyst = T,
                  mulname = "all-qc",
                  ...)
    
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-k",default = 8,type= "integer", help = "聚类数")
  parser$add_argument("-s",default = 3,type= "integer", help = "稀疏参数")
  parser$add_argument("-r",default = 1,type= "integer", help = "距离半径")
  parser$add_argument("-mk","--mulk",default = 15,type= "integer", help = "多样本聚类数")
  parser$add_argument("-ms","--muls",default = 3,type= "integer", help = "多样本稀疏参数")
  parser$add_argument("-mr","--mulr",default = 1,type= "integer", help = "多样本距离半径")
  parser$add_argument("-al","--allanalyst",default = "T",help = "是否对all进行sscc")

  args <- parser$parse_args()

  writeinfo()
  createdir(filename = "./sample/cluster",linkdir = T)
  args$allanalyst <- as.logical(args$allanalyst)
  
  result <- do.call(what = MulanaCluster,args = args)

  writeinfo(endtime = T)
}
