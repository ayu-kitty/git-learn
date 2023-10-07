#!/opt/conda/bin/Rscript

#' 对目录下imzml进行信息统计
#'
#' @param slidename 选择玻片名
#' @param imzmlpath 数据路径
#' @param moderange 正负离子模式
#' @param ... 见[getimzmlbasicinfo()]
#'
#' @export
imzmlinfo <- function(slidename = NULL,
                      imzmlpath = "./sample/imzml/",
                      moderange = c("neg", "pos"),
                      savepath = imzmlpath,
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
      
      mse <- readdata(mse)
      
      getimzmlbasicinfo(mse = mse,
                        mode = mode,
                        filename = paste0(savepath,"/",samplename[i], "-", mode,".xlsx"),
                        ...)
      getimzmlbasicinfo(mse = mse[,1:10],
                        mode = mode,
                        filename = paste0(savepath,"/",samplename[i], "-", mode,"-bg.xlsx"),
                        ...)
      
      gc(reset = TRUE)
      
    }
  }
}

#' 获取imzml文件基础信息
#'
#' @param mse imzml读取的mse数据
#' @param freq.min 最低频率,滤波范围(出现比例少于此参数的峰会被剔除)
#' @param tolerance 峰校正参数,检测粉于基准对齐的公差,默认20(有用传参)
#' @param units 分辨率单位
#' @param rm.zero 布尔,是否删除平均强度为零的特征
#' @param filename 保存文件名
#' 
#' @export
getimzmlbasicinfo <- function(mse,
                              freq.min = 0.1,
                              tolerance = 20,
                              units = "ppm",
                              rm.zero = T,
                              filename = "test.xlsx",
                              mode = "neg",
                              waters = F,
                              ...){
  
  setCardinalBPPARAM(MulticoreParam(workers = 3))
  
  mse <- readdata(mse)
  
  if(dim(mse)[2] > 1000){
    print("~数据大于1000个像素点，减少数据至1000像素点进行运算")
    # 随机选点1000个像素点
    selectnum <- sample(x = 1:dim(mse)[2],size = 1000)
    # 排序
    selectnum <- selectnum[order(selectnum)]
    # 取出新的mse
    mse2 <- mse[,selectnum]
  }else{
    mse2 <- mse
  }
  
  print(mse2)
  
  source(file = packagepath(path = "script/spatial/peakAlign.R"))
  
  #峰对齐,峰筛选,cardinal包的peakalign函数有问题,所以新写了一个peakalign2,为了waters扩大了窗口为3
  mse3 <- mse2 %>%
    peakAlign2(tolerance = tolerance, units = units) %>%
    peakFilter(freq.min = freq.min, rm.zero = rm.zero) %>%
    process()
  
  mse3 <- mse2 %>%
    peakBin(ref = mz(mse3),type= "height",tolerance = tolerance, units = units) %>% 
    process()
  
  # 筛选基线
  freqintensity <- 100
  # 筛选基线以上的表达值
  intensityselect <- featureApply(mse3, calzeroratio,intensity = freqintensity) >= freq.min
  # 提取基线以上的表达值
  mse3 <- subsetFeatures(mse3,intensityselect)
  
  # 取表达量均值
  samplemeandata <- featureApply(mse3, mean)
  
  # 取表达量极大值
  samplemaxdata <- featureApply(mse3, outliermax)
  
  # 保存infodata
  infodata <- data.frame(mz = mz(mse3),
                         samplemean = samplemeandata,
                         samplemaxdata = samplemaxdata)
  
  # 计算高于基线强度的值的比例
  infodata[,"zeroratiodata"] <- featureApply(mse3, calzeroratio,freqintensity)
  # 计算rsd
  infodata[,"rsd"] <- featureApply(mse3, function(x){sd(x)/mean(x)})
  infodata2 <- data.frame("信息" = c("最大强度","平均强度","mz数","mz"),
                          "数据" = c(max(samplemeandata),mean(samplemeandata),sum(intensityselect),NA))
  
  if(mode == "neg"){
    mz <- c(89.02446,114.93512,115.89073,118.03994,121.0284,
            123.90097,129.05465,130.93011,134.86409,138.56808,
            139.89592,141.01955,143.10674,148.9233,159.88605,
            180.03038,183.00389,189.05682,239.05913,255.23281,
            277.14428,283.26372,554.2615,617.25701,639.23851)
  }else{
    mz <- c(140.06842,171.09941,173.07873,187.09429,201.0736,
            207.09948,261.13093,273.16745,279.15858,301.14034,
            317.11503,437.19315,453.16686,556.2771,579.2921)
  }
  
  if(length(mz) > 0){
    
    if(waters){
      if(mode == "neg"){
        mz2 <- 554.26202254
      }else{
        mz2 <- 556.27657454
      }
      mse_sample <- mzAlign3(object = mse2, ref = mz2, refmz = mz,
                             tolerance = tolerance, units = units)
    }else{
      mse_sample <- mse2 %>%
        peakBin(ref = mz,type= "height",tolerance = tolerance, units = units) %>% 
        process()
    }
    
    infodata3 <- data.frame("信息" = mz(mse_sample),
                            "数据" = featureApply(mse_sample, mean))
    
    infodata2 <- rbind(infodata2,infodata3)
  }
  
  savexlsx1(data = infodata2,filename = filename,sheet = "基本信息")
  savexlsx1(data = infodata,filename = filename,sheet = "离子信息")
  
  setCardinalBPPARAM(SerialParam())
}
