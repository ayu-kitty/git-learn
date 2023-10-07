#!/opt/conda/bin/Rscript

#' 对目录下imzml进行空代数据归一化
#'
#' @param slidename 选择玻片名
#' @param imzmlpath 数据路径
#' @param moderange 正负离子模式
#' @param ... 见[getimzmlNormalization()]
#'
#' @export
imzmlNormalization <- function(slidename = NULL,
                               imzmlpath = "./sample/imzml/",
                               moderange = c("neg", "pos"),
                               savepath = "./sample/imzml-pre/",
                               ...) {
  suppressMessages(library("Cardinal"))
  
  # 正负离子循环
  for (mode in moderange) {
    
    # 读数据列表
    if (is.null(slidename)) {
      filename <- list.files(path = imzmlpath,
                             pattern = paste0("-", mode, ".imzML$"),
                             full.names = T,
                             recursive = T)
    } else {
      filename <- paste0(imzmlpath,"/",slidename, "-", mode, ".imzML")
      filename <- filename[file.exists(filename)]
    }
    
    # 判断数据是否存在
    if (length(filename) == 0) {
      print(paste0("在", mode, "模式下",imzmlpath,"中未找到imzml文件"))
      next
    }
    
    # 样本名处理
    samplename <- gsub(pattern = paste0("-", mode, ".imzML"),
                       replacement = "",
                       x = basename(filename))
    
    referenceintensity <- NA
      
    for ( i in 1:length(samplename)) {
      
      # 归一化计算
      referenceintensity <- getimzmlNormalization(filename = filename[i],
                                                  samplename = samplename[i],
                                                  mode = mode,
                                                  savepath = savepath,
                                                  referenceintensity = referenceintensity,
                                                  ...)

      gc(reset = TRUE)
      
    }
  }
}


#' 对单一的imzml进行空代数据归一化
#' @param filename 文件名
#' @param samplename 样本名
#' @param mode 正负离子模式
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param referencemz 归一化模式
#' @param tolerance waters矫正参数,此处默认10
#' @param savepeakpath peak保存路径,默认"./sample/peak/"
#' @param savepath imzml预处理保存路径,默认"./sample/imzml-pre/"
#' @param normethod 归一化方式,默认none(还有tic和rms)
#' @param referencemethod 归一化模式,默认"pixel"
#' @param norreferencemzneg reference归一化模式下neg选择的离子 89.02319,#277.14428
#' @param norreferencemzpos reference归一化模式下pos选择的离子 301.14034,#317.1140
#' @param referenceintensity reference归一化表达强度,默认NA
#' @param stripdeal 是否不进行条纹处理,默认F,开启条纹处理
#' @param imagetype 图像保存格式
#' @param asp 比例,默认1
#' @param waters 是否为waters仪器,默认F
#' 
#' @export
getimzmlNormalization <- function(filename = NULL,
                                  samplename = NULL,
                                  mode = "neg",
                                  # 读取参数
                                  mass.range = NULL,
                                  resolution = 5,
                                  units = "ppm",
                                  # 峰对齐参数
                                  referencemz = T,
                                  tolerance = 10,
                                  savepeakpath = "./sample/peak/",
                                  savepath = "./sample/imzml-pre/",
                                  # 归一化方式
                                  normethod = "none",
                                  referencemethod = "pixel",
                                  norreferencemzneg = 89.02319,#277.14428,
                                  norreferencemzpos = 301.14034,#317.1140,
                                  referenceintensity = NA,
                                  # 条纹处理
                                  stripdeal = F,
                                  imagetype = "jpg",
                                  asp = 1,
                                  waters = F,
                                  number = 100000,
                                  ...){
  suppressMessages(library("Cardinal"))
  setCardinalBPPARAM(MulticoreParam(workers = 3))
  
  print(paste0(filename,"处理中"))
  
  # if(waters){
  #   mse <- readimzml(filename = filename,
  #                    mass.range = mass.range,
  #                    resolution = 0.002,
  #                    units = "mz",
  #                    attach.only = F,
  #                    ...)
  # }else{
  # 读数据
    mse <- readimzml(filename = filename,
                     mass.range = mass.range,
                     resolution = resolution,
                     units = units,
                     attach.only = F,
                     ...) 
  # }

  print(mse)
  print(paste0("~使用分辨率为:",resolution(mse),names(resolution(mse))))
  
  # 峰对齐
  # 获取referencemz-raw的数据
  if (referencemz & file.exists(paste0(savepeakpath,"referencemz-raw-", mode, ".rds"))) {
    print("使用已获取mz值")
    mz <- readdata(filename = paste0(savepeakpath,"referencemz-raw-", mode, ".rds"))
    
    # mse <- mse %>%
    #   peakAlign(tolerance = tolerance, units = units, ref = mz) %>%
    #   process()
    # mse <- pull(mse, as.matrix = TRUE)
    
    if(!is.na(tolerance)){
      if(tolerance <= 0){
        tolerance <- NA
      }
    }
    
    # waters的归一矫正mz改变为554和556
    if(waters){
      # if(mode == "neg"){
      #   mz2 <- negrefmz
      #   norreferencemzneg <- mz2
      # }else{
      #   mz2 <- posrefmz
      #   norreferencemzpos <- mz2
      # }
      
      mz2 <- readdata(filename = paste0(savepeakpath,"referencemz-raw-waters-", mode, ".rds"))
      
      # source(file = packagepath(path = "script/spatial/mzAlign.R"))
      # 峰对齐
      mse <- mzAlign3(object = mse, ref = mz2, refmz = mz,
                      tolerance = tolerance, units = units,number = number,
                      savepath = savepath)
      
      # mse <- mse %>%
      #   mzAlign2(ref = mz2,tolerance = 0.1,units ="mz") %>%
      #   peakBin(ref = mz,tolerance = tolerance, units = units) %>%
      #   # smoothSignal(method="gaussian") %>%
      #   process()
      
      # if(mode == "neg"){
      #   spectra(mse)[,] <- t(featureApply(mse,function(x){x[x < 20] <- 0;x})) 
      # }else{
      #   spectra(mse)[,] <- t(featureApply(mse,function(x){x[x < 40] <- 0;x})) 
      # }
      
      # 图片保存
      # 保存内标成像图,
      imzmlimage(filename = mse,
                 savepath = savepath,
                 mapname = paste0(samplename, "-image-", mode),
                 mapmz = T,
                 mz = mz2,
                 imagetype = imagetype,
                 asp = asp)
      # 保存内标质谱图
      imzmlplot(filename = mse,
                mapmz = T,
                mz = mz2,
                savepath = savepath,
                mapname = paste0(samplename,"-plot-",mode),
                imagetype = imagetype)
      
    }else{
      # HF 峰对齐
      mse <- mse %>%
        peakBin(ref = mz,type= "height",tolerance = tolerance, units = units) %>% 
        process()
    }

    print(mse)
    print(paste0("~使用峰对齐参数为:",tolerance))
  }else{
    stop("无对齐mz信息")
  }
  
  # 内标归一化
  if(normethod == "reference"){
    # 内标mz获取
    norreferencemz <- ifelse(mode=="neg",norreferencemzneg,norreferencemzpos)
    # 所有mz获取
    allmz <- mz(mse)
    num <- NULL
    for ( k in 1:length(norreferencemz)) {
      # 获取聚类归一化内标最近的mz
      mzlist <- abs(allmz-norreferencemz[k])
      # 最近mz的位数
      num2 <- which.min(mzlist)
      # 最近mz的前一位和后一位
      num2 <- c(num2-1,num2,num2+1)
      # 获取mz列表
      num <- c(num,num2)
    }
    # mz列表排序
    num <- num[order(num)]
    # mz 列表去重
    num <- num[!duplicated(num)]
    # 获取列表对应的mz
    mzrange <- allmz[num]
    # 获取列表对应mz的mse
    mse_sample2 <- mse[num,]
    
    # mz列表取均值最大的mz
    diffmz <- mzrange[which.max(featureApply(mse_sample2,mean))]
    # allmz与diffmz的距离绝对值
    mzlist <- abs(allmz-diffmz)
    # 最近的mz
    num <- which.min(mzlist)
    # mz与归一化的norreferencemz的最小距离(绝对值)
    normzdelta <-min(abs(norreferencemz-allmz[num]))
    
    # 计算归一内标mz理论值和实际值的ppm
    normzppm <- normzdelta/allmz[num]*1000000
    # 导入num位的光谱文件 
    refintensitydata <- spectra(mse)[num,]
    # 计算normz的rsd
    normzsd <- sd(refintensitydata)
    # normz的均值
    normzintensity <- mean(refintensitydata)
    # normz的中位数
    normzmedianintensity <- median(refintensitydata)
    # normz的rsd
    normzrsd <- normzsd/normzintensity
    # 求和表达量>10/像素点总数
    normzratio <- sum(refintensitydata > 10)/dim(mse)[2]
    # 打印
    print(paste0("~normz:",allmz[num]))
    print(paste0("~normzppm:",normzppm))
    print(paste0("~normzintensity:",normzintensity))
    print(paste0("~normzmedianintensity:",normzmedianintensity))
    print(paste0("~normzratio:",normzratio))
    print(paste0("~normzrsd:",normzrsd))
    
    # 保存内标的质谱图的jpg
    imzmlimage(filename = mse,
               savepath = savepath,
               mapname = paste0(samplename, "-ref-", mode),
               mapmz = T,
               mz = allmz[num],
               imagetype = "jpg")
    
    # 内标表达强度低于200/2000,覆盖度低于 0.1,视为无内标
    if(mode=="neg"){
      if(normzintensity < 200 | normzratio < 0.1){
        stop("~未找到背景离子")
      }
    }else{
      if(normzintensity < 2000 | normzratio < 0.1){
        stop("~未找到背景离子")
      }
    }
    
    # 如果未传入内标强度,则使用normz表达强度
    if(is.na(referenceintensity)){
      referenceintensity <- normzmedianintensity
    }
    
    
    if(referencemethod == "sample"){
      # 以sample样本做归一
      mse <- Cardinal::pull(mse, as.matrix = TRUE)
      spectradata <- as.matrix(spectra(mse))
      spectradata <- spectradata/normzmedianintensity*referenceintensity
      spectra(mse) <- spectradata
    }else{
      # 以reference为标准归一
      mse <- mse %>%
        normalize(method = "reference",feature = num,scale = referenceintensity) %>%
        process()
    }
    
    # 去除内标mz的mse
    mse <- mse[-num,]
  }else if(normethod == "tic")
    # TIC归一化
    {
    # 如果未传入内标强度,则求和总强归一
    if(is.na(referenceintensity)){
      referenceintensity <- median(pixelApply(mse,sum))
    }
    
    # tic归一
    mse <- mse %>%
      Cardinal::normalize(method = "tic",tic = referenceintensity) %>%
      process()
  }else if(normethod == "rms"){
    
    # 如果未传入内标强度,则求和总强中位数归一
    if(is.na(referenceintensity)){
      referenceintensity <- median(pixelApply(mse,sum))
    }
    
    mse <- mse %>%
      Cardinal::normalize(method = "rms",rms = referenceintensity) %>%
      process()
  }else{
    print("~不进行归一化")
  }
  
  # 条纹处理
  if(stripdeal){
    mse <- imzmlstripdeal(mse = mse)
  }
  
  # 数据保存
  filename <- paste0(savepath,"/",samplename, "-", mode,".imzML")
  saveimzml(data = mse,
            filename = filename)
  
  # 图片保存
  imzmlimage(filename = filename,
             savepath = savepath,
             mapname = paste0(samplename, "-", mode),
             imagetype = imagetype,
             asp = asp)
  
  setCardinalBPPARAM(SerialParam())
  
  return(referenceintensity)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-i","--imzmlpath",default = "./sample/imzml/", help = "imzml原始文件路径,默认./sample/imzml/")
  parser$add_argument("-s","--savepath",default = "./sample/imzml-pre/", help = "imzml数据保存路径,默认./sample/imzml-pre/")
  parser$add_argument("-sp","--savepeakpath",default = "./sample/peak/", help = "mz数据保存路径,默认./sample/peak/")
  parser$add_argument("-m","--massrange",default =NULL,type= "double",help = "mz范围,默认NULL",
                      nargs = 2,dest = "mass.range")
  parser$add_argument("-re","--resolution",default = 5, type= "double",help = "分辨率,默认5")
  parser$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  parser$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  parser$add_argument("-a","--asp",default = 1, type= "double",help = "长宽分辨率比")
  
  parser$add_argument("-sl","--slidename",default = NULL,help = "选择处理玻片",nargs = "+")
  
  # 峰对齐及筛选参数
  parser$add_argument("-tr","--tolerance",default = 10, type= "double",help = "峰对齐参数,默认5")
  parser$add_argument("-w","--waters",default = F,help = "是否使用waters仪器",action ='store_true')
  parser$add_argument("-nb","--number",default = 100000, type= "integer",help = "分段矫正参数")
  
  # 内标归一化
  parser$add_argument("-nm","--normethod",default = "none", help = "是否不进行内标归一化处理,不进行填写none",
                      choices = c("none","reference","tic"))
  parser$add_argument("-nrn","--norreferencemzneg",default = 89.02319, type= "double", help = "reference归一化模式下,选择的负离子",nargs = "+")
  parser$add_argument("-nrp","--norreferencemzpos",default = 301.14034, type= "double", help = "reference归一化模式下,选择的正离子",nargs = "+")
  
  # 条纹处理
  parser$add_argument("-sd","--stripdeal",default = F, help = "是否不进行条纹处理",action='store_true',
                      dest = "stripdeal")
  
  args <- parser$parse_args()
  
  writeinfo()
  createdir(filename = args$savepath,linkdir = T)
  
  result <- do.call(what = imzmlNormalization,args = args)
  
  writeinfo(endtime = T)
}
