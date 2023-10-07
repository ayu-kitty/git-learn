#!/opt/conda/bin/Rscript

#' 对imzml数据批量压缩聚类
#'
#' @export
AllCompressAndCluster <- function(imzmlpath = "./sample/imzml/",
                                  moderange = c("neg", "pos"),
                                  imzmlmoudle = CompressAndCluster,
                                  samplename = NULL,
                                  ...) {
  
  print("进行压缩聚类,初步区分背景样本区域")
  # 初步压缩聚类(函数位于3.1.2-5-1 imzmlanacall.R)
  imzmlanacall(imzmlpath = imzmlpath,
               moderange = moderange,
               imzmlmoudle = CompressAndCluster,
               samplename = samplename,
               ...)
}


#' 对imzml数据压缩聚类
#'
#' @param samplename 样本名称
#' @param mode 正负离子模式
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param r 距离半径
#' @param k 聚类数量
#' @param s 稀疏参数
#' @param iter.max 聚类循环次数
#' @param imzmlpath 文件路径
#' @param savepath 保存路径
#' @param imagetype 图像格式
#' @param rebg 逻辑，是否剔除背景离子
#' @param bgimzmlpath 背景数据路径
#' @param bgintensity 背景离子强度
#' @param bgfreq 背景离子表达比例
#' @param ... 见[clusterlimage()]
#'
#' @export
CompressAndCluster <- function(samplename = "A",
                               mode = "neg",
                               mass.range = NULL,
                               resolution = 100,
                               units = "ppm",
                               r = 1,
                               k = 10,
                               s = 9,
                               iter.max = 30,
                               imzmlpath = "./sample/imzml/",
                               savepath = "./sample/compress/",
                               imagetype = "jpg",
                               rebg = F,
                               bgimzmlpath = "./sample/compress-raw/",
                               bgintensity = 100,
                               bgfreq = 0.05,
                               multipleintensity = 2,
                               attach.only = F,
                               asp = 1,
                               ...) {
  suppressMessages(library("Cardinal"))
  setCardinalBPPARAM(MulticoreParam(workers = 3))
  
  # print(paste0(imzmlpath, samplename, "-", mode, ".imzML", "文件---处理中"))
  # read
  mse <- readimzml(filename = paste0(imzmlpath, samplename, "-", mode, ".imzML"),
                   attach.only = attach.only,
                   mass.range = mass.range, 
                   resolution = resolution, 
                   units = units)
  
  print("~~数据压缩中~~")
  
  #  峰过滤,检测粉于基准对齐的公差为100,压缩数据????????为什么是100
  mse <- mse %>%
    peakAlign(tolerance = 100, units = "ppm") %>%
    peakFilter(freq.min = 0.1, rm.zero = T) %>%
    # normalize(method = "tic") %>%
    process()
  
  mse_proc <- mse
  mse_proc <- pull(mse_proc, as.matrix = TRUE)
  # spectra(mse_proc) <- as.matrix(spectra(mse_proc))
  
  if(rebg){
    if(file.exists(paste0(bgimzmlpath, "bg_data-",mode, ".imzML"))){
      print("~~背景离子剔除中~~")
      # 读取bgdata 的imzml 
      bgmse <- readMSIData(paste0(bgimzmlpath, "bg_data-",mode, ".imzML"),
                           attach.only = attach.only,
                           mass.range = mass.range, resolution = resolution, units = units)
      
      if(dim(bgmse)[2] > 10000){
        print("~背景数据大于10000个像素点，减少数据至10000像素点进行运算")
        # 随机选点1000个
        selectnum <- sample(x = 1:dim(bgmse)[2],size = 10000)
        selectnum <- selectnum[order(selectnum)]
        bgmse <- bgmse[,selectnum]
      }
      
     # bgdata峰过滤 ,公差100
      bgmse <- bgmse %>%
        peakAlign(tolerance = 100, units = "ppm",ref = mz(mse_proc)) %>%
        process()
      
      # 用cardinal包实现bgdata 提取
      bgratio <- featureApply(bgmse, calzeroratio, bgintensity)
      bgmeanintensity <- featureApply(bgmse, mean)
      allmeanintensity <- featureApply(mse_proc, mean)
      select <- (bgratio < bgfreq) | (allmeanintensity > bgmeanintensity*multipleintensity)
      
      while(sum(select) < 100){
        print("~~剔除背景离子后小于100个离子~~")
        bgfreq <- bgfreq*1.1
        print(paste0("~~重新设置背景频率为",bgfreq,"~~"))
        select <- (bgratio < bgfreq) | (allmeanintensity > bgmeanintensity*multipleintensity)
      }
      
      mse_proc <- mse_proc[select,]
      
      print(paste0("~~",sum(select),"个mz被保留~~"))
      
    }else{
      print(paste0("~~",bgimzmlpath, mode, ".imzML文件不存在~~"))
    }
  }
  
  setCardinalBPPARAM(SerialParam())
  print("~~数据压缩完成~~")
  
  # print("~~压缩数据保存中~~")
  filename <- paste0(savepath,"/",samplename,"-",mode,".imzML")
  saveimzml(data = mse_proc,
            filename = filename)
  
  imzmlimage(filename = filename,
             savepath = savepath,
             mapname = paste0(samplename, "-", mode),
             imagetype = imagetype,
             asp = asp)

  # 聚类，背景删除
  print("~~数据聚类中~~")
  ssc <- spatialShrunkenCentroids(mse_proc, 
                                  r = r, 
                                  k = k, 
                                  s = s, 
                                  iter.max = iter.max)
  
  # print("聚类数据保存中")
  saverds(data = list(sscc = ssc),
          filename = paste0(savepath,"/",samplename, "-", mode, ".rds"))
  
  # print("聚类结果绘图中")
  imzmlclusterimage(filename = paste0(savepath,"/",samplename, "-", mode, ".rds"),
                    savepath = savepath,
                    mapname = paste0(samplename, "-cluster-", mode),
                    imagetype = imagetype,
                    asp = asp)
  
  return("~~处理完成~~")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-i","--imzmlpath",default = "./sample/imzml/", help = "imzml原始文件路径,默认./sample/imzml/")
  parser$add_argument("-s","--savepath",default = "./sample/compress/", help = "imzml数据保存路径,默认./sample/compress/")
  parser$add_argument("-sa","--samplename",default = NULL, help = "选择样本进行分析",nargs = "+")
  parser$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  parser$add_argument("-a","--asp",default = 1, type= "double",help = "长宽分辨率比")
  parser$add_argument("-m","--massrange",default = NULL,type= "double",help = "mz范围,默认NULL",
                      nargs = 2,dest = "mass.range")
  parser$add_argument("-re","--resolution",default = 100, type= "double",help = "分辨率,默认100")
  parser$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  
  parser$add_argument("-rb","--rebg",default = F, help = "是否剔除背景离子",action='store_true')
  parser$add_argument("-b","--bgimzmlpath",default = "./sample/compress-raw/", help = "背景imzml数据保存路径,默认./sample/compress-raw/")
  parser$add_argument("-bi","--bgintensity",default = 100,  type= "double",help = "背景离子识别强度")
  parser$add_argument("-bf","--bgfreq",default = 0.05,  type= "double",help = "背景离子表达范围")
  parser$add_argument("-mi","--multipleintensity",default = 2,  type= "double",help = "样本离子比背景离子表达倍数")
  args <- parser$parse_args()
  
  writeinfo()
  
  mulargs <- do.call(what = AllCompressAndCluster,args = args)
  
  writeinfo(endtime = T)
  
}
