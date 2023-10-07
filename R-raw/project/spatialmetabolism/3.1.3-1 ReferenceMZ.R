#!/opt/conda/bin/Rscript

#' 对目录下imzml获取定量mz
#'
#' @param sample 选择样本名
#' @param areapath 选区路径
#' @param moderange 正负离子模式
#' @param maxnum 最大样本数
#' @param seed 随机选择种子
#' @param ... 见[GetimzmlReferenceMZ()]
#' @param imzmlpath imzml数据路径
#' @param delsample 不运行的样本名
#'
#' @export
imzmlReferenceMZ <- function(sample = NULL,
                             imzmlpath = "./sample/imzml/",
                             areapath = "./sample/select/",
                             moderange = c("neg", "pos"),
                             maxnum = 10,
                             seed = 1111,
                             delsample = c("bg_data","qc_data"),
                             ...) {
  
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
    # 样本名
    samplename <- gsub(pattern = ".*\\/", replacement = "", x = samplename2)
    
    # 样本名拆分
    slidename <- gsub(pattern = "\\/.*", replacement = "", x = samplename2)
    
    # 玻片名(imzml数据名称)
    realslidename <- list.files(path = imzmlpath,pattern = paste0("-",mode,".imzML$"))
    realslidename <- gsub(pattern = paste0("-",mode,".imzML"),
                          replacement = "",
                          x = realslidename)
    
    # 获的需要处理的扣完bg的样本和bgdata名
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
    
    # sample数量判断
    if(!is.null(sample)){
      samplename2 <- samplename2[samplename %in% sample]
      samplename <- samplename[samplename %in% sample]
    }else if(length(samplename2) > maxnum){
      # 大于maxnum(10)数量的样本,随机选择样本进行后续处理
      set.seed(seed =  seed)
      selectnum <- sample(1:length(samplename2),size = maxnum)
      samplename2 <- samplename2[selectnum]
      samplename <- samplename[selectnum]
    }
    
    if(length(samplename2) > 0){
      
      slidename <- gsub(pattern = "\\/.*", replacement = "", x = samplename2)
      filename <- paste0(imzmlpath,"/",slidename,"-",mode,".imzML")
      areards <- paste0(areapath,"/",mode,"/",samplename2,".rds")
      
      # 单个mz获取定量信息(见函数)
      GetimzmlReferenceMZ(filename = filename,
                          mode = mode,
                          areards = areards,
                          ...)
      
    }else{
      print(paste0("在",areapath,"目录下未找到", mode, "模式的样本文件"))
    }
    
  }
}

#' 对单个imzml获取定量mz
#'
#' @param filename 文件名 
#' @param areards 区域rds路径
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param tolerance 对齐参数
#' @param freq.min 0值筛选,滤波范围(出现比例少于此参数的峰会被剔除)
#' @param rm.zero 去除0值
#' @param mode 正负离子模式
#' @param savepeakpath peak存储路径
#' @param number 每number个像素点拆分进行分段矫正,默认3000
#' @param ...
#'
#' @export
GetimzmlReferenceMZ <- function(filename = NULL,
                                areards = NULL,
                                mode = "neg",
                                # 读取参数
                                mass.range = NULL,
                                resolution = 5,
                                units = "ppm",
                                # 峰对齐参数
                                tolerance = 10,
                                freq.min = 0.01,
                                rm.zero = T,
                                # 峰筛选参数
                                savepeakpath = "./sample/peak/",
                                asp = 1,
                                waters = F,
                                number = 100000,
                                # negrefmz = c(554.26202254),
                                # posrefmz = c(556.27657454),
                                negrefmz = c(89.0244,124.0074,133.0142,171.1391,255.233,327.233,554.262,714.5079,834.5283),
                                posrefmz = c(104.1075,147.1128,149.0233,156.042,239.1642,301.141,556.2766,578.2585,606.2942,734.5694,772.5238,798.5401,826.5705,844.5334),
                                ...) {
  suppressMessages(library("Cardinal"))
  setCardinalBPPARAM(MulticoreParam(workers = 3))
  
  if(!is.na(tolerance)){
    if(tolerance <= 0){
      tolerance <- NA
    }
  }
  
  mse_sample <- readimzml(filename = filename,
                          mass.range = mass.range,
                          resolution = resolution,
                          units = units,
                          area = T,
                          areards = areards,
                          addy = T,
                          ...)
  # mse_sample <<- mse_sample
  print(mse_sample)
  print(paste0("~使用分辨率为:",resolution(mse_sample),names(resolution(mse_sample))))
  print(paste0("~使用峰对齐参数为:",tolerance,names(resolution(mse_sample))))
  
  if(dim(mse_sample)[2] > 50000){
    print("~数据大于50000个像素点，减少数据至50000像素点进行运算")
    selectnum <- sample(x = 1:dim(mse_sample)[2],size = 50000)
    selectnum <- selectnum[order(selectnum)]
    mse_sample <- mse_sample[,selectnum]
  }
  
  # water文件的mz 读取list
  watersmzfile <- list.files(path = "./sample/waterspeak/",pattern = mode,full.names = T)
  if(waters & length(watersmzfile)>0){
    
    mz <- readdata(watersmzfile[1])
    mz <- mz[,1]
    mz <- mz[order(mz)]
    
    # 内标mz提取
    if(mode == "neg"){
      mz2 <- negrefmz
    }else{
      mz2 <- posrefmz
    }
    
    # source(file = packagepath(path = "script/spatial/mzAlign.R"))
    # 质量轴矫正,ref = mz2(可以考虑mz2参数外置),函数见3.1.3-2-2-mzAlign
    
    mse_sample <- mulmzAlign3(object = mse_sample, ref = mz2, refmz = mz,
                              tolerance = tolerance, units = units,number = number,
                              savepath = savepeakpath)
    mse_sample <- mse_sample %>%
      peakFilter(freq.min = freq.min, rm.zero = rm.zero) %>%
      process()
    
    # mse_sample <- mse_sample %>%
    #   mzAlign2(ref = mz2,tolerance = 0.1,units ="mz") %>%
    #   peakAlign(ref = mz,tolerance = tolerance, units = units) %>%
    #   peakFilter(freq.min = freq.min, rm.zero = rm.zero) %>%
    #   process()
    
  }else if(waters & length(watersmzfile)==0){
    
    # 内标mz提取
    if(mode == "neg"){
      mz2 <- negrefmz
    }else{
      mz2 <- posrefmz
    }
    
    # mse_sample <<- mse_sample
    
    mse_sample <- mulmzAlign3(object = mse_sample, ref = mz2, refmz = NULL,
                              tolerance = tolerance, units = units,number = number,
                              savepath = savepeakpath)
    
    source(file = packagepath(path = "script/spatial/peakAlign.R"))
    mse_sample <- mse_sample %>%
      peakAlign2(tolerance = tolerance, units = units) %>%
      peakFilter(freq.min = freq.min, rm.zero = rm.zero) %>%
      process()
    
  }else{
    
    source(file = packagepath(path = "script/spatial/peakAlign.R"))
    # 正常的hf sample提取(对齐参考)
    mse_sample <- mse_sample %>%
      peakAlign2(tolerance = tolerance, units = units) %>%
      peakFilter(freq.min = freq.min, rm.zero = rm.zero) %>%
      process()
    
  }
  
  if(mode == "neg"){
    freqintensity <- 50
  }else{
    freqintensity <- 100
  }
  
  # 大于初筛基线强度的mz数量占总比的freq.min(0.01)(TF)
  intensityselect <- featureApply(mse_sample, calzeroratio,intensity = freqintensity) >= freq.min
  
  # 取出mse符合上一条标准的的mse
  mse_sample <- mse_sample[intensityselect,]
  
  # 取出mz
  mz <- mz(mse_sample)
  
  # mz个数
  num <- length(mz)
  
  print(paste0("~",mode,"模式下离子数量为：",num))
  
  # 保存对齐后的mz数据
  saverds(data = mz,
          filename = paste0(savepeakpath,"referencemz-raw-", mode, ".rds"))
  
  if(waters){
    
    # intensity <- featureApply(mse_sample, mean)
    # mz <- mz(mse_sample)[intensity > 2000]
    # for ( mz3 in mz2) {
    #   mz <- mz[(mz < (mz3 - 0.1)) | (mz > (mz3 + 0.1))]
    # }
    # mz <- c(mz,mz2)
    # mz2 <- mz[order(mz)]
    saverds(data = mz2,
            filename = paste0(savepeakpath,"referencemz-raw-waters-", mode, ".rds"))
    
  }
  
  setCardinalBPPARAM(SerialParam())
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-i","--imzmlpath",default = "./sample/imzml/", help = "imzml原始文件路径,默认./sample/imzml/")
  parser$add_argument("-a","--areapath",default = "./sample/select/", help = "选区文件路径,默认./sample/select/")
  parser$add_argument("-sp","--savepeakpath",default = "./sample/peak/", help = "mz数据保存路径,默认./sample/peak/")
  parser$add_argument("-m","--massrange",default = NULL,type= "double",help = "mz范围,默认NULL",
                      nargs = 2,dest = "mass.range")
  parser$add_argument("-re","--resolution",default = 5, type= "double",help = "分辨率,默认5")
  parser$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  parser$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  
  parser$add_argument("-sa","--sample",default = NULL,help = "选择对齐样本")
  parser$add_argument("-ds","--delsample",default = c("bg_data","qc_data"),help = "删除对齐样本",nargs = "+")
  # 峰对齐及筛选参数
  parser$add_argument("-tr","--tolerance",default = 10, type= "double",help = "峰对齐参数,默认5")
  parser$add_argument("-fm","--freqmin",default = 0.01, type= "double",help = "最小表达范围比例,默认0.01",
                      dest = "freq.min")
  
  parser$add_argument("-mn","--maxnum",default = 10, type= "integer",help = "最大样本数")
  parser$add_argument("-se","--seed",default = 1111, type= "integer",help = "样本选择随机种子")
  parser$add_argument("-w","--waters",default = F,help = "是否使用waters仪器",action ='store_true')
  parser$add_argument("-nb","--number",default = 100000,type= "integer",help = "分段矫正参数")
  parser$add_argument("-nrf","--negrefmz",default = c(89.0244,124.0074,133.0142,171.1391,255.233,327.233,554.262,714.5079,834.5283), type= "double",help = "负离子的质量轴校正",nargs = "+")
  parser$add_argument("-prf","--posrefmz",default = c(104.1075,147.1128,149.0233,156.042,239.1642,301.141,556.2766,578.2585,606.2942,734.5694,772.5238,798.5401,826.5705,844.5334), type= "double",help = "正离子的质量轴校正",nargs = "+")
  
  args <- parser$parse_args()
  
  writeinfo()
  createdir(filename = args$savepeakpath,linkdir = T)
  
  result <- do.call(what = imzmlReferenceMZ,args = args)
  
  writeinfo(endtime = T)
  
}
