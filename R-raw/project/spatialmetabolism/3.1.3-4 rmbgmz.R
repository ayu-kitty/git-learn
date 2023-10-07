#!/opt/conda/bin/Rscript

#' 剔除背景及低表达离子
#'
#' @param sample 选择样本名
#' @param areapath 选区路径
#' @param moderange 正负离子模式
#' @param maxnum 最大样本数
#' @param seed 随机选择种子
#' @param imzmlpath imzml数据路径
#' @param delsample 不运行的样本名
#' @param ... 见[Getimzmlrmbgmz()]
#'
#' @export
imzmlrmbgmz <- function(sample = NULL,
                        imzmlpath = "./sample/imzml-pre/",
                        areapath = "./sample/select/",
                        moderange = c("neg", "pos"),
                        maxnum = 10,
                        seed = 1111,
                        delsample = c("bg_data","qc_data"),
                        ...) {
  suppressMessages(library("Cardinal"))
  
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
    
    # 玻片名
    slidename <- gsub(pattern = "\\/.*", replacement = "", x = samplename2)
    realslidename <- list.files(path = imzmlpath,pattern = paste0("-",mode,".imzML$"))
    realslidename <- gsub(pattern = paste0("-",mode,".imzML"),
                          replacement = "",
                          x = realslidename)
    # 样本名确认
    samplename2 <- samplename2[slidename %in% realslidename]
    samplename <- samplename[slidename %in% realslidename]
    
    # bg名获取
    bgsamplename <- samplename2[samplename %in% "bg_data"]
    
    # 超过maxnum (10)的选择10个进行背景离子扣除
    if(length(bgsamplename) > maxnum){
      print(paste0("背景超过",maxnum,"个,将进行随机选择"))
      set.seed(seed =  seed)
      selectnum <- sample(1:length(bgsamplename),size = maxnum)
      bgsamplename <- bgsamplename[selectnum]
    }
    print("选择以下背景:")
    print(bgsamplename)
    
    # 从sample列表中删除背景选区及QC
    if(!is.null(delsample)){
      for(delsample2 in delsample) {
        samplename2 <- samplename2[!grepl(pattern = paste0("^",delsample2,".*"),x = samplename)]
        samplename <- samplename[!grepl(pattern = paste0("^",delsample2,".*"),x = samplename)]
      }
    }
    
    
    if(!is.null(sample)){
      # 确认sample是否在samplename列表中,是否有对应数据
      samplename2 <- samplename2[samplename %in% sample]
      samplename <- samplename[samplename %in% sample]
    }else if(length(samplename2) > maxnum){
      # 超过maxnum(10)个样本的sample也随机选maxnum(10)个
      print(paste0("样本超过",maxnum,"个,将进行随机选择"))
      set.seed(seed =  seed)
      selectnum <- sample(1:length(samplename2),size = maxnum)
      samplename2 <- samplename2[selectnum]
      samplename <- samplename[selectnum]
    }
    
    # 样本名重复报错
    if (any(duplicated(samplename))) {
      print(samplename)
      stop("样本名有重复，请核查")
    }
    
    # 打印筛选和确认过的sample列表
    print("选择以下样本:")
    print(samplename2)
    
    
    if(length(samplename2) > 0){
      
      slidename <- gsub(pattern = "\\/.*", replacement = "", x = samplename2)
      filename <- paste0(imzmlpath,"/",slidename,"-",mode,".imzML")
      areards <- paste0(areapath,"/",mode,"/",samplename2,".rds")
      
      # bg区域数据提取
      if(length(bgsamplename) > 0){
        bgslidename <- gsub(pattern = "\\/.*", replacement = "", x = bgsamplename)
        bgfilename <- paste0(imzmlpath,"/",bgslidename,"-",mode,".imzML")
        bgareards <- paste0(areapath,"/",mode,"/",bgsamplename,".rds")
      }else{
        bgfilename <- NULL
        bgareards <- NULL
      }
      
      #对单个imzml剔除背景及低表达离子
      Getimzmlrmbgmz(filename = filename,
                     areards = areards,
                     bgfilename = bgfilename,
                     bgareards = bgareards,
                     mode = mode,
                     ...)
      
    }else{
      print(paste0("在",areapath,"目录下未找到", mode, "模式的样本文件"))
    }
    
  }
}


#' 对单个imzml剔除背景及低表达离子
#'
#' @param filename 文件名路径
#' @param areards 文件对应rsd路径
#' @param mode 正负离子模式
#' @param mass.range mass 范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param attach.only readmzml参数,详见readMSIData()
#' @param savepeakpath peak数据存储路径,默认"./sample/peak/"
#' @param negmaxintensityratio neg离子最大表达强度1000
#' @param negmeanintensityratio neg离子平均表达强度200
#' @param negfreqintensityratio neg离子范围表达强度200
#' @param posmaxintensityratio pos离子最大表达强度1000
#' @param posmeanintensityratio pos离子平均表达强度200
#' @param posfreqintensityratio pos离子范围表达强度200
#' @param freq 离子表达范围,默认0.01
#' @param minsmoothnum 每个样本最小连续表达数,默认16
#' @param maxsmoothnum 1个样本时,最大连续表达数,默认100
#' @param rebg 是否剔除背景区域,默认T
#' @param bgfilename 背景文件路径
#' @param bgareards 背景文件rsd路径
#' @param negbgintensityratio neg背景离子表达强度
#' @param posbgintensityratio pos背景离子表达强度
#' @param bgratio 背景离子表达范围默认0.2
#' @param bgtomeansampleratio 样本与背景离子的均值比,默认1.2
#' @param bgtomaxsampleratio 样本与背景离子最大值比,默认1.5
#' @param # qualitative 是否定性
#' @param # saveannopath 定性结果保存路径,默认"./sample/qualitative/" 
#' @param # mapmz 是否绘制mz成像图
#' @param # asp 图像比例,默认1
#' @param rmouter 边缘处理,默认F,旧流程pos外轮廓删除
#' 
#' @export
#' 
Getimzmlrmbgmz <- function(filename = NULL,
                           areards = NULL,
                           mode = "neg",
                           # 读取参数
                           mass.range = NULL,
                           resolution = 5,
                           units = "ppm",
                           attach.only = F,
                           savepeakpath = "./sample/peak/",
                           # 根据强度筛选
                           negmaxintensityratio = 1000,
                           negmeanintensityratio = 100,
                           negfreqintensityratio = 100,
                           posmaxintensityratio = 1000,
                           posmeanintensityratio = 100,
                           posfreqintensityratio = 100,
                           freq = 0.01,
                           # 根据连续表达数
                           minsmoothnum = 16,
                           maxsmoothnum = 100,
                           # 根据背景筛选
                           rebg = T,
                           bgfilename = NULL,
                           bgareards = NULL,
                           negbgintensityratio = 100,
                           posbgintensityratio = 100,
                           bgratio = 0.2,
                           bgtomeansampleratio = 1.2,
                           bgtomaxsampleratio = 1.5,
                           # 根据定性结果筛选
                           # qualitative = T,
                           # saveannopath = "./sample/qualitative/",
                           # mapmz = F,
                           # asp = 1,
                           rmouter = F,
                           ...) {
  suppressMessages(library("Cardinal"))
  setCardinalBPPARAM(MulticoreParam(workers = 3))
  
  # mse读取(所有sample),多个run
  mse_sample <- readimzml(filename = filename,
                          mass.range = mass.range,
                          resolution = resolution,
                          units = units,
                          area = T,
                          areards = areards,
                          addy = T,
                          attach.only = attach.only,
                          ...)
  
  # 异常值去除后的最大表达强度
  samplemaxintensity <- max(featureApply(mse_sample, outliermax))
  
  # 正负离子最大&平均&范围表达强度传参 & 背景离子表达强度传参
  if(mode == "neg"){
    maxintensityratio  <- negmaxintensityratio 
    meanintensityratio <- negmeanintensityratio
    freqintensityratio <- negfreqintensityratio
    bgintensityratio <- negbgintensityratio
  }else if(mode == "pos"){
    maxintensityratio  <- posmaxintensityratio 
    meanintensityratio <- posmeanintensityratio
    freqintensityratio <- posfreqintensityratio
    bgintensityratio <- posbgintensityratio
  }else{
    stop(paste0("无",mode,"模式"))
  }
  
  # 如果表达强度传参是小数,则 表达强度=所有样本最大表达强度*表达强度传参， 如果不是小数则直接赋值
  if(maxintensityratio <= 1){maxintensity <- maxintensityratio*samplemaxintensity
  }else{maxintensity <- maxintensityratio}
  if(meanintensityratio <= 1){meanintensity <- meanintensityratio*samplemaxintensity
  }else{meanintensity <- meanintensityratio}
  if(freqintensityratio <= 1){freqintensity <- freqintensityratio*samplemaxintensity
  }else{freqintensity <- freqintensityratio}
  if(bgintensityratio <= 1){bgintensity <- bgintensityratio*samplemaxintensity
  }else{bgintensity <- bgintensityratio}

  print(mse_sample)
  
  print("计算样本均值")
  # 计算异常值去除后的均值
  samplemeandata <- featureApply(mse_sample, outliermean)
  # samplemeandata <- featureApply(mse_sample, mean)
  # samplemeandata <- featureApply(mse_sample, median)
  print("计算样本最大值")
  # 计算样本表达强度最大值
  samplemaxdata <- featureApply(mse_sample, outliermax)
  
  # info信息提取,此处mz就是未筛选的mz数
  infodata <- data.frame(mz = mz(mse_sample),
                         sampleselect = T,
                         samplemeandata = samplemeandata,
                         samplemaxdata = samplemaxdata)
  
  print(paste0("初始有",sum(infodata[,"sampleselect"]),"个mz"))
  
  # 根据背景强度剔除mz
  if(rebg){
    if(length(bgfilename) > 0){
      print("背景剔除流程")
      # bg区域mse数据读取
      bgmse <- readimzml(filename = bgfilename,
                         mass.range = mass.range,
                         resolution = resolution,
                         units = units,
                         area = T,
                         areards = bgareards,
                         ...)
      
      if(dim(bgmse)[2] > 50000){
        print("~数据大于50000个像素点，减少数据至50000像素点进行运算")
        # 随机筛选50000像素点处理
        bgmse <- bgmse[,sample(x = 1:dim(bgmse)[2],size = 50000)]
      }
      
      # 计算bg区域的mean和max
      print("背景强度计算")
      infodata[,"bgmeandata"] <- featureApply(bgmse, mean)
      infodata[,"bgmaxdata"] <- featureApply(bgmse, outliermax)
      
      print("背景强度比例计算")
      # 背景强度的覆盖度计算(比例,bg区域占样本总区域的比例) calzeror atio:  sum(x > intensity)/length(x)
      infodata[,"bgratiodata"] <- featureApply(bgmse,calzeroratio,bgintensity)
      
      # 背景强度覆盖度>=背景离子表达范围 为T
      infodata[,"bgratioselect"] <- (infodata$bgratiodata >= bgratio)
      # 表达强度均值>=背景离子均值*样本与背景离子的均值比 为T
      infodata[,"sampleintensitymeanselect"] <- (infodata$samplemeandata >= infodata$bgmeandata*bgtomeansampleratio)
      #表达强度最大值>=背景离子最大值*样本与背景离子最大值比 为T
      infodata[,"sampleintensitymaxselect"] <- (infodata$samplemaxdata >= infodata$bgmaxdata*bgtomaxsampleratio)
      
      # select mz标注(非bgratioselect区域,sampleintensitymeanselect 和 sampleintensitymaxselect)
      infodata[,"sampleselectinbg"] <- !infodata$bgratioselect
      infodata[,"sampleselectinbg"] <- infodata$sampleselectinbg | infodata$sampleintensitymeanselect | infodata$sampleintensitymaxselect
      # sampleselect 为去除背景离子后的mz
      infodata[,"sampleselect"] <- infodata$sampleselectinbg & infodata$sampleselect
  
      print(paste0("~经过背景离子扣除后,有",sum(infodata[,"sampleselect"]),"个mz被保留"))

    }else{
      print("~未找到背景文件，不进行背景强度筛选")
    }
  
  }
  
  # 根据样本强度剔除mz
  print("样本强度及分布剔除流程")
  # 获取每个像素点的x,y坐标,run信息(样本来源)作为group
  xdata <- coord(mse_sample)$x
  ydata <- coord(mse_sample)$y
  group <- as.character(run(mse_sample))
  
  # 样本强度比例计算
  print("样本强度比例计算")
  # 大于均值表达强度传参的 为T
  infodata[,"selectmean"] <- (infodata$samplemeandata >= meanintensity)
  # 大于最大值表达强度传参 为T
  infodata[,"selectmax"] <- (infodata$samplemaxdata >= maxintensity)
  # 范围表达强度占比计算
  infodata[,"zeroratiodata"] <- featureApply(mse_sample, calzeroratio,freqintensity)
  # 选择强度占比大于离子表达范围的mz
  infodata[,"zeroratioselect"] <- (infodata$zeroratiodata >= freq)
  # 满足selectmean,selectmax,zeroratioselect 任意一个为T的 做保留
  infodata[,"intensityselect"] <- infodata[,"selectmean"] | infodata[,"selectmax"] | infodata[,"zeroratioselect"]
  
  print("样本连续表达数计算")
  # if(mode=="pos"){
  #   infodata[,"smoothnum"] <- featureApply(mse_sample,smoothintensitynum,x= xdata,y = ydata,
  #                                          minintensity = freqintensity,rmouter = rmouter)
  # 
  # }else{
  #   infodata[,"smoothnum"] <- featureApply(mse_sample,smoothintensitynum,x= xdata,y = ydata,
  #                                          minintensity = freqintensity,rmouter = F)
  # }
  # 
  # infodata[,"numselect"] <-  (infodata[,"smoothnum"] >= minsmoothnum)
  
  # infodata[,"numselect"] <- featureApply(mse_sample,smoothintensitynum,x= xdata,y = ydata,
  #                                        minintensity = freqintensity,rmouter = rmouter,minsmoothnum = minsmoothnum)
  
  # 连续表达筛选详见smoothintensitynum和mulsmoothintensitynum函数(下方)
  #smoothintensitynum：单样本表达强度筛选判定     mulsmoothintensitynum： 所有样本连续表达强度筛选判定
  # 对于单个样本,mz在连续表达超过maxsmoothnum(100)个像素点中表达的视为连续表达
  infodata[,"numselect"] <- featureApply(mse_sample,smoothintensitynum,x= xdata,y = ydata,
                                         minintensity = freqintensity,rmouter = rmouter,minsmoothnum = maxsmoothnum)|
    # 对于所有样本,mz在连续表达超过minsmoothnum(16)个像素点中表达的视为连续表达
    featureApply(mse_sample,mulsmoothintensitynum,x= xdata,y = ydata,group = group,
                 minintensity = freqintensity,rmouter = rmouter,minsmoothnum = minsmoothnum)
  
  # 最终选择,根据上述的bgsampel进行的sampleselect,根据表达强度进行的intensityselect,根据连续表达进行的numselect筛选的交集,作为最后的选择
  infodata[,"sampleselect"] <- infodata[,"sampleselect"] & infodata[,"intensityselect"]  & infodata[,"numselect"]
  
  # 根据定性结果筛选
  # if(qualitative & file.exists(paste0(saveannopath,"/adducts-",mode,".rds"))){
  #   qualitativedata <- readRDS(file = paste0(saveannopath,"/adducts-",mode,".rds"))
  #   qualitativedata <- qualitativedata[,c("mz","adducts","formula","isotope","ppm")]
  #   qualitativedata2 <- qualitativedata[!(qualitativedata$mz %in% infodata$mz),]
  #   qualitativedata <- qualitativedata[!(qualitativedata$formula %in% qualitativedata2$formula),]
  #   infodata2 <- merge(infodata,qualitativedata,by = "mz",all.x = T)
  #   
  #   if(mapmz){
  #     mzlist <- infodata2[,c("mz","sampleselect","formula")]
  #     mzlist <- mzlist[!is.na(mzlist$formula),]
  #     mzlist <- mzlist[mzlist$sampleselect,]
  #     mzlist2 <- list()
  #     for (formula in unique(mzlist$formula)) {
  #       mzlist2 <- c(mzlist2,list(mzlist[mzlist$formula == formula,"mz"]))
  #     }
  #     imzmlimage(filename = filename[1],
  #                mass.range = mass.range,
  #                resolution = resolution,
  #                units = units,
  #                area = T,
  #                areards = areards[1],
  #                savepath = paste0(saveannopath,"/map/",mode,"/"),
  #                mapname = paste0("sample-", mode),
  #                imagetype = "jpg",
  #                mapmz = T,
  #                mz = mzlist2,
  #                normalize.image = "linear",
  #                asp = asp,
  #                ...)
  #   }
  #   
  #   infodata2[,"qualitative"] <- ((!is.na(infodata2$adducts)) & is.na(infodata2$isotope))
  #   infodata2[,"sampleselect"] <- infodata2[,"sampleselect"]&infodata2[,"qualitative"]
  #   infodata3 <- infodata2[infodata2$sampleselect,]
  #   infodata3 <- infodata3[order(infodata3$samplemeandata,decreasing = T),]
  #   infodata3 <- infodata3[!duplicated(infodata3$formula),]
  #   infodata2[!(infodata2$mz %in% infodata3$mz),"sampleselect"] <- F
  #   infodata <- infodata2
  # }
  
  # 根据最终筛选的mz,获得mz列表
  mz <- mz(mse_sample)[infodata[,"sampleselect"]]
  
  # 保存为referencemz.rds
  saverds(data = mz,
          filename = paste0(savepeakpath,"referencemz-", mode, ".rds"))
  
  # 保存筛选依据的infodata
  savexlsx1(data = infodata,
            filename = paste0(savepeakpath,"样本离子统计.xlsx"),
            sheet = mode)
  
  print(paste0("~经过筛选后,有",sum(infodata[,"sampleselect"]),"个mz被保留"))
  
  setCardinalBPPARAM(SerialParam())
}

#' 单样本连续表达强度筛选判定
#'
#' @param data mse_sample
#' @param x x坐标
#' @param y y坐标
#' @param group 样本分组(sample名)
#' 
#' @export
mulsmoothintensitynum <-function(data,x,y,group,...){
  
  allresult <- c()
  
  # 单样拆分
  for (sample in unique(group)) {
    data2 <- data[group == sample]
    x2 <- x[group == sample]
    y2 <- y[group == sample]  
    # 单样进行表达强度筛选判定
    result <- smoothintensitynum(data = data2,
                                 x = x2,
                                 y = y2,
                                 ...)
    allresult <- c(allresult,result)
  }
  
  # 返回符合要求的连续表达强度大于50%的视为连续表达
  return(sum(allresult)/length(allresult) > 0.5)
  
}

#' 所有样本表达强度筛选判定
#'
#' @param data mse_sample
#' @param x x坐标
#' @param y y坐标
#' @param minintensity 离子范围表达强度(从negfreqintensityratio和posfreqintensityratio及后续计算传入)
#' @param rmouter 是否边缘处理
#' @param minsmoothnum 连续表达标准(minsmoothnum or maxsmoothnum 16 or 100)
#' 
#' @export
smoothintensitynum <- function(data,x,y,minintensity = 0,rmouter = F,minsmoothnum = 9){

  # data <- spectra(mse)[1,]
  alldata <- data.frame(intensity = data,x = x,y = y)

  if(any(data!=0)){
    alldata[alldata$intensity < quantile(data[data!=0],0.25),"intensity"] <- 0 
  }else{
    return(F)
  }
  
  # 旧版pos边缘扣除,常规为F,不调用
  if(rmouter){
    if(max(x)-min(x) > 50){
      for ( i in unique(y)) {
        alldata[alldata$y == i & alldata$x == max(alldata[alldata$y == i,"x"]),"intensity"] <- 0
        alldata[alldata$y == i & alldata$x == (max(alldata[alldata$y == i,"x"])-1),"intensity"] <- 0
        alldata[alldata$y == i & alldata$x == min(alldata[alldata$y == i,"x"]),"intensity"] <- 0
        alldata[alldata$y == i & alldata$x == (min(alldata[alldata$y == i,"x"])+1),"intensity"] <- 0
      }
    }else{
      for ( i in unique(y)) {
        alldata[alldata$y == i & alldata$x == max(alldata[alldata$y == i,"x"]),"intensity"] <- 0
        alldata[alldata$y == i & alldata$x == min(alldata[alldata$y == i,"x"]),"intensity"] <- 0
      }
    }
    
    if(max(y)-min(y) > 50){
      for ( i in unique(x)) {
        alldata[alldata$x == i & alldata$y == max(alldata[alldata$x == i,"y"]),"intensity"] <- 0
        alldata[alldata$x == i & alldata$y == (max(alldata[alldata$x == i,"y"])-1),"intensity"] <- 0
        alldata[alldata$x == i & alldata$y == min(alldata[alldata$x == i,"y"]),"intensity"] <- 0
        alldata[alldata$x == i & alldata$y == (min(alldata[alldata$x == i,"y"])+1),"intensity"] <- 0
      }
    }else{
      for ( i in unique(x)) {
        alldata[alldata$x == i & alldata$y == max(alldata[alldata$x == i,"y"]),"intensity"] <- 0
        alldata[alldata$x == i & alldata$y == min(alldata[alldata$x == i,"y"]),"intensity"] <- 0
      }
    }
  }
  # 小于等于离子范围表达强度的置为0,大于离子表达强度的置为1
  alldata[alldata$intensity <= minintensity,"intensity"] <- 0
  alldata[alldata$intensity > minintensity,"intensity"] <- 1
  # 取出大于1的部分
  alldata <- alldata[alldata$intensity == 1,]
  # alldata <<- alldata
  
  # 大于1为空的返回F
  if(dim(alldata)[1] == 0){
    return(F)
  }
  
  # 创建一个x长,y宽的matrix,作为样本的坐标映射
  alldata2 <- matrix(NA, 
                     nrow = max(y), ncol = max(x),
                     dimnames = list(1:max(y),1:max(x)))
  
  # 判断如果满足离子表达强度筛选的为1
  for ( i in 1:dim(alldata2)[1]) {
    alldata2[i,as.character(alldata[alldata$y == i,"x"])] <- 1
  }
  # 空值置为0
  alldata2[is.na(alldata2)] <- 0
  
  arealist <- list()
  maxArea <- 0
  
  # 判断相邻点
  for ( k in 1:dim(alldata)[1]) {
    # 循环满足离子强度范围筛选的所有强度值,获取x,y坐标
    i <- alldata[k,"y"]
    j <- alldata[k,"x"]
    if (alldata2[i,j] == 1) {
      alldata2[i,j] <- 0
      areanum <- 1
      # 获取上下左右的4个点的坐标
      arealist <- c(arealist,list(c(i-1,j),c(i+1,j),c(i,j-1),c(i,j+1)))
      while (length(arealist) > 1) {
        i0 <- arealist[[1]][1]
        j0 <- arealist[[1]][2]
        # 如果周边一个点为1,那么areanum+1
        if((i0 <= 0) | (j0 <= 0) | (i0 > dim(alldata2)[1] ) | (j0 > dim(alldata2)[2])){
        }else if(alldata2[i0,j0] == 1){
          alldata2[i0,j0] <- 0
          areanum <- areanum + 1
          if(areanum >= minsmoothnum){
            return(T)
          }
          arealist <- c(arealist,list(c(i0-1,j0),c(i0+1,j0),c(i0,j0-1),c(i0,j0+1)))
        }
        
        arealist[[1]] <- NULL
      }
      if(areanum > maxArea){
        maxArea <- areanum
      }
    }
  }
  # 如果最大连续值大于minsmoothnum(100 or 16),则返回T,否则F
  if(maxArea >= minsmoothnum){
    return(T)
  }else{
    return(F)
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-i","--imzmlpath",default = "./sample/imzml-pre/", help = "imzml原始文件路径,默认./sample/imzml-pre/")
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
  
  parser$add_argument("-mn","--maxnum",default = 10, type= "integer",help = "最大样本数")
  parser$add_argument("-se","--seed",default = 1111, type= "integer",help = "样本选择随机种子")
  
  # 背景筛选参数
  parser$add_argument("-nair","--negmaxintensityratio",default = 1000, type= "double",help = "离子最大表达强度")
  parser$add_argument("-nmir","--negmeanintensityratio",default = 100, type= "double",help = "离子平均表达强度")
  parser$add_argument("-nnfi","--negfreqintensityratio",default = 100, type= "double",help = "离子范围表达强度")
  parser$add_argument("-pair","--posmaxintensityratio",default = 1000, type= "double",help = "离子最大表达强度")
  parser$add_argument("-pmir","--posmeanintensityratio",default = 100, type= "double",help = "离子平均表达强度")
  parser$add_argument("-pnfi","--posfreqintensityratio",default = 100, type= "double",help = "离子范围表达强度")
  parser$add_argument("-fq","--freq",default = 0.01, type= "double",help = "离子表达范围")
  parser$add_argument("-nbir","--negbgintensityratio",default = 100, type= "double",help = "离子背景表达强度")
  parser$add_argument("-pbir","--posbgintensityratio",default = 100, type= "double",help = "离子背景表达强度")
  parser$add_argument("-br","--bgratio",default = 0.2, type= "double",help = "离子背景表达范围")
  parser$add_argument("-bm","--bgtomeansampleratio",default = 1.2,type= "double",help = "离子背景与样本均值比")
  parser$add_argument("-ba","--bgtomaxsampleratio",default = 1.5,type= "double",help = "离子背景与样本最大值比")
  parser$add_argument("-mnn","--minsmoothnum",default = 16, type= "integer",help = "每个样本最小连续表达数,默认9")
  parser$add_argument("-mxn","--maxsmoothnum",default = 100, type= "integer",help = "仅有一个样本最大连续表达数,默认100")
  parser$add_argument("-sap","--saveannopath",default = "./sample/qualitative/", help = "定性保存路径,默认./sample/qualitative/")
  parser$add_argument("-mm","--mapmz",default = F, help = "是否绘制相关成像图",action='store_true')
  
  args <- parser$parse_args()
  
  writeinfo()
  
  createdir(filename = args$savepeakpath,linkdir = T)
  result <- do.call(what = imzmlrmbgmz,args = args)
  
  writeinfo(endtime = T)
}
