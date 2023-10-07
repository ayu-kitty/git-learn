#!/opt/conda/bin/Rscript

#' 空代批量样本预处理
#'
#' @param sample 选择样本名
#' @param areapath 选区路径
#' @param moderange 正负离子模式
#' @param imzmlpath imzml数据路径
#' @param delsample 不运行的样本名
#' @param savepeakpath peak信息路径
#' @param saveannopath 定性信息路径
#' @param ... 见[PreDeal()]
#' @param qualitative 是否生成定性结果
#' @param deldup 是否删除重复物质
#' @param referencemzpath 不生成定性结果时获取mz的路径
#'
#' @export
MulPreDeal <- function(sample = NULL,
                       imzmlpath = "./sample/imzml-pre/",
                       areapath = "./sample/select/",
                       moderange = c("neg", "pos"),
                       delsample = c("bg_data"),
                       qualitative = T,
                       savepeakpath = "./sample/peak/",
                       saveannopath = "./sample/qualitative/",
                       deldup = T,
                       referencemzpath = "./sample/peak/",
                       referencemzname = "referencemz",
                       ...) {
  if(qualitative){
    qualitativedata <- readdata(filename = paste0(saveannopath,"/adducts.xlsx"))
    if("sampleselect" %in% colnames(qualitativedata)){
      peakdata <- getfinalqualitative_2023(moderange = moderange,
                                           saveannopath =  saveannopath,
                                           deldup = deldup)
    }else{
      peakdata <- getfinalqualitative(moderange = moderange,
                                      savepeakpath =  savepeakpath,
                                      saveannopath =  saveannopath,
                                      deldup = deldup)
    }
  }else{
    peakdata <- NULL
    for (mode in moderange) {
      if(file.exists(paste0(referencemzpath,"/",referencemzname,"-",mode,".rds"))){
        peakdata2 <- list(readdata(paste0(referencemzpath,"/",referencemzname,"-",mode,".rds")))
        names(peakdata2) <- mode
        peakdata <- c(peakdata,peakdata2)
      }
    }
    # print(peakdata)
  }
  
  # 正负离子循环
  for (mode in moderange) {
    
    filename <- list.files(path = paste0(areapath, mode),
                           pattern = ".rds$",
                           full.names = F,
                           recursive = T)
    # 获取样品名
    samplename2 <- gsub(pattern = ".rds",
                        replacement = "",
                        x = filename)
    samplename <- gsub(pattern = ".*\\/", replacement = "", x = samplename2)
    
    slidename <- gsub(pattern = "\\/.*", replacement = "", x = samplename2)
    realslidename <- list.files(path = imzmlpath,pattern = paste0("-",mode,".imzML$"))
    realslidename <- gsub(pattern = paste0("-",mode,".imzML"),
                          replacement = "",
                          x = realslidename)
    samplename2 <- samplename2[slidename %in% realslidename]
    samplename <- samplename[slidename %in% realslidename]
    
    # 删除背景选区及QC
    if(!is.null(delsample)){
      for(delsample2 in delsample) {
        samplename2 <- samplename2[!grepl(pattern = paste0("^",delsample2,".*"),x = samplename)]
        samplename <- samplename[!grepl(pattern = paste0("^",delsample2,".*"),x = samplename)]
      }
    }
    
    if (any(duplicated(samplename))) {
      stop("样本名有重复，请核查")
    }
    
    if(!is.null(sample)){
      samplename2 <- samplename2[samplename %in% sample]
      samplename <- samplename[samplename %in% sample]
    }
    
    if(length(samplename2) > 0){
      
      for (i in seq_len(length(samplename2))) {

        PreDeal(slidename = strsplit(x = samplename2[i], split = "\\/")[[1]][1],
                samplename = strsplit(x = samplename2[i], split = "\\/")[[1]][2],
                mode = mode,
                areapath = areapath,
                imzmlpath = imzmlpath,
                peakdata = peakdata[[mode]],
                ...)

      }
      
    }else{
      print(paste0("在",areapath,"目录下未找到", mode, "模式的样本文件"))
    }
    
  }
}


#' 空代单样本预处理
#'
#' @param slidename 玻片名称
#' @param samplename 样本名称
#' @param mode 正负离子模式
#' @param imzmlpath 数据存储位置
#' @param areapath 区域存储位置
#' @param rotatepath 旋转存储位置 
#' @param peakpath peak存储位置
#' @param presavepath 预处理结果存储路径
#' @param finalsavepath 最终结果存储路径
#' @param ...
#'
#' @export
PreDeal <- function(slidename = NULL,
                    samplename = NULL,
                    mode = "neg",
                    peakdata = NULL,
                    imzmlpath = "./sample/imzml-pre/",
                    areapath = "./sample/select/",
                    rotatepath = paste0(areapath,"/rotate/"),
                    presavepath = "./sample/pre/",
                    finalsavepath = "./sample/final/",
                    mvcoord = T,
                    stripdeal = T,
                    imagetype = "jpg",
                    asp = 1,
                    ...) {
  
  print(paste0(samplename, "-", mode, "运行中"))
  suppressMessages(library("Cardinal"))
  mse <- readimzml(filename = paste0(imzmlpath, slidename, "-", mode, ".imzML"),
                   attach.only = F,
                   ...)
  mz <- mz(mse)
  mz <- format(mz , nsmall = 5, trim = T)
  peakdata <- format(peakdata , nsmall = 5, trim = T)
  base::print(peakdata[!(peakdata %in% mz)])
  peakdata <- peakdata[peakdata %in% mz]
  base::print(paste0("mz数量为:",length(peakdata)))
  mse <- subsetFeatures(mse,mz = as.numeric(peakdata))
  
  print("读取样本区域信息")
  area <- readrds(filename = paste0(areapath, mode, "/", slidename, "/", samplename, ".rds"))
  
  if(file.exists(paste0(rotatepath, mode, "/", slidename, "/", samplename, ".rds"))){
    print("读取样本旋转信息")
    rotate <- readrds(filename = paste0(rotatepath, mode, "/", slidename, "/", samplename, ".rds"))
  }else{
    print("未读取到样本旋转信息,使用默认值")
    rotate <- list(mirrorx = F,
                   mirrory = F,
                   trans = F)
  }
  
  print("空代数据处理中")
  mse$sample <- area
  
  xmax <- max(coord(mse)$x[area]) + 5
  xmin <- min(coord(mse)$x[area]) - 5
  ymax <- max(coord(mse)$y[area]) + 5
  ymin <- min(coord(mse)$y[area]) - 5
  
  mse <- subsetPixels(mse,
                      x <= xmax & x >= xmin,
                      y <= ymax & y >= ymin)
  
  if(mvcoord){
    coord(mse)$x <- coord(mse)$x - (min(coord(mse)$x) - 1)
    coord(mse)$y <- coord(mse)$y - (min(coord(mse)$y) - 1)
  }

  if(rotate$mirrorx){
    coord(mse)$x <- max(coord(mse)$x) + min(coord(mse)$x) - coord(mse)$x
  }
  if(rotate$mirrory){
    coord(mse)$y <- max(coord(mse)$y) + min(coord(mse)$y) - coord(mse)$y
  }
  if(rotate$trans){
    coordxy <- coord(mse)$x
    coord(mse)$x <- coord(mse)$y
    coord(mse)$y <- coordxy
  }
  
  print("处理背景区域")
  mse <-  Cardinal::pull(mse, as.matrix = TRUE)
  spectra(mse)[,!mse$sample] <- 0
  
  spectra(mse)[,] <- t(featureApply(mse,outliertomax))

  # 条纹处理
  if(stripdeal){
    mse <- imzmlstripdeal(mse = mse,saverawdata = F)
  }
  
  spectra(mse)[,!mse$sample] <- 0
  
  print("数据保存中")
  
  mse2 <- mse
  mse2 <- imzmlfill(mse2)
  
  mse$samplename <- samplename
  mse$slidename <- slidename
  run(mse) <- samplename
  
  # 保存给客户结果
  featureinfo <- data.frame(sample = mse2$sample,
                            x = coord(mse2)$x,
                            y = coord(mse2)$y,
                            run = run(mse2),
                            samplename = samplename,
                            slidename = slidename,
                            mode = mode)
  
  filename <- paste0(presavepath,"/",samplename, "-", mode,".imzML")
  
  saveimzml(data = mse2,
            filename = filename)
  saverds(data = featureinfo,
          filename = paste0(presavepath,"/",samplename, "-", mode, ".rds"))
  
  imzmlimage(filename = filename,
             savepath = presavepath,
             mapname = paste0(samplename, "-", mode),
             imagetype = imagetype,
             asp = asp)
  
  # 保存最终结果
  mse <- subsetPixels(mse, sample)
  featureinfo <- data.frame(sample = mse$sample,
                            x = coord(mse)$x,
                            y = coord(mse)$y,
                            run = run(mse),
                            samplename = samplename,
                            slidename = slidename,
                            mode = mode)
  
  filename <- paste0(finalsavepath,"/",samplename, "-", mode,".imzML")
  
  saveimzml(data = mse,
            filename = paste0(finalsavepath,"/",samplename, "-", mode))
  saverds(data = featureinfo,
          filename = paste0(finalsavepath,"/",samplename, "-", mode, ".rds"))
  
  imzmlimage(filename = filename,
             savepath = finalsavepath,
             mapname = paste0(samplename, "-", mode),
             imagetype = imagetype,
             asp = asp)
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-i","--imzmlpath",default = "./sample/imzml-pre/", help = "imzml原始文件路径,默认./sample/imzml-pre/")
  parser$add_argument("-ap","--areapath",default = "./sample/select/", help = "选区文件路径,默认./sample/select/")
  parser$add_argument("-sp","--savepeakpath",default = "./sample/peak/", help = "mz数据保存路径,默认./sample/peak/")
  parser$add_argument("-m","--massrange",default = NULL,type= "double",help = "mz范围,默认NULL",
                      nargs = 2,dest = "mass.range")
  parser$add_argument("-re","--resolution",default = 5, type= "double",help = "分辨率,默认5")
  parser$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  parser$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  parser$add_argument("-a","--asp",default = 1, type= "double",help = "长宽分辨率比")
  
  parser$add_argument("-sa","--sample",default = NULL,help = "选择对齐样本")
  parser$add_argument("-ds","--delsample",default = c("bg_data"),help = "删除对齐样本",nargs = "+")
  parser$add_argument("-sap","--saveannopath",default = "./sample/qualitative/", help = "定性结果保存路径,默认./sample/qualitative/")
  
  # 结果保存
  parser$add_argument("-psp","--presavepath",default = "./sample/pre/", help = "给客户数据保存路径,默认./sample/pre/")
  parser$add_argument("-fsp","--finalsavepath",default = "./sample/final/", help = "最终数据保存路径,默认./sample/final/")
  parser$add_argument("-mc","--mvcoord",default = T, help = "是否不进行坐标调整",action='store_false')
  # 条纹处理
  parser$add_argument("-nsd","--nostripdeal",default = T, help = "是否不进行条纹处理",action='store_false',
                      dest = "stripdeal")
  
  # 定性结果
  parser$add_argument("-nq","--noqualitative",default = T, help = "是否不用定性结果做峰提取",action='store_false',
                      dest = "qualitative")
  parser$add_argument("-rn","--referencemzname",default = "referencemz", help = "不用定性结果后使用的峰对齐数据")
  
  args <- parser$parse_args()
  
  writeinfo()
  createdir(filename = args$presavepath,linkdir = T)
  createdir(filename = args$finalsavepath,linkdir = T)
  
  result <- do.call(what = MulPreDeal,args = args)
  
  writeinfo(endtime = T)

}
