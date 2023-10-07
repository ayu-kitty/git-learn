#!/opt/conda/bin/Rscript

#' 背景数据提取
#'
#' @param imzmlpath 数据路径
#' @param clusterfrom 背景聚类来源
#' @param infopath 登记单路径
#' @param moderange 正负离子模式
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param savepath 保存路径
#' @param ... [见clusterlimage()]
#'
#' @export
imzmlbgdata <- function(imzmlpath = "./sample/imzml/",
                        clusterfrom = "./sample/compress-raw/",
                        infopath = "项目登记单.xlsx",
                        moderange = c("neg","pos"),
                        savepath = clusterfrom,
                        areapath = "./sample/select/",
                        mass.range = NULL,
                        resolution = 5,
                        units = "ppm",
                        attach.only = F,
                        imagetype = "jpg",
                        asp = 1,
                        ...){
  suppressMessages(library("Cardinal"))
  
  if (!file.exists(infopath)) {
    stop("不存在项目登记单.xlsx")
  }
  
  # 项目信息获取
  bginfo <- readdata(filename = infopath, sheet = "背景信息")
  
  for (mode in moderange) {
    bginfomode <- bginfo[bginfo$模式 == mode,]
    
    if(dim(bginfomode)[1] > 0){
      
      j <- 1
      mse <- readimzml(filename = paste0(imzmlpath, "/", bginfomode[j,"玻片名"], "-", mode, ".imzML"),
                       attach.only = attach.only,
                       mass.range = mass.range, 
                       resolution = resolution, 
                       units = units)
      
      ssc <- readrds(filename = paste0(clusterfrom, "/", bginfomode[j,"玻片名"], "-", mode, ".rds"))
      
      ssccluster <- ssc$sscc$class[[1]]
      arearange <- as.character(bginfomode[j,"聚类"])
      arearange <- as.numeric(strsplit(x = arearange, split = "\\+")[[1]])
      arearange2 <- ssccluster %in% arearange
      
      savearea(slidename = bginfomode[j,"玻片名"],
               samplename = "bg_data",
               mode = mode,
               area = arearange2,
               clusterpath = clusterfrom,
               savepath = areapath,
               asp = asp)
      
      mse <- mse[,arearange2]
      msebg <- mse
      
      if(dim(bginfomode)[1] > 1){
        for ( j in 2:dim(bginfomode)[1]) {
          
          mse <- readimzml(filename = paste0(imzmlpath, "/", bginfomode[j,"玻片名"], "-", mode, ".imzML"),
                           attach.only = attach.only,
                           mass.range = mass.range, 
                           resolution = resolution, 
                           units = units)
          
          ssc <- readrds(filename = paste0(clusterfrom, "/", bginfomode[j,"玻片名"],"-", mode, ".rds"))
          
          ssccluster <- ssc$sscc$class[[1]]
          arearange <- as.character(bginfomode[j,"聚类"])
          arearange <- as.numeric(strsplit(x = arearange, split = "\\+")[[1]])
          arearange2 <- ssccluster %in% arearange
          
          savearea(slidename = bginfomode[j,"玻片名"],
                   samplename = "bg_data",
                   mode = mode,
                   area = arearange2,
                   clusterpath = clusterfrom,
                   savepath = areapath,
                   asp = asp)
          
          mse <- mse[,arearange2]
          msebg <- BiocGenerics::cbind(msebg, mse)
          
        }
      }
      
      run(msebg) <- mode
      
      if(dim(msebg)[2] > 10000){
        print("~背景数据大于10000个像素点，减少数据至10000像素点进行运算")
        selectnum <- sample(x = 1:dim(msebg)[2],size = 10000)
        selectnum <- selectnum[order(selectnum)]
        msebg <- msebg[,selectnum]
      }
      
      saveimzml(data = msebg,
                filename = paste0(savepath,"/","bg_data-",mode,".imzML"))
      
    }else{
      print(paste0(mode,"模式下未选择背景"))
    }
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-i","--imzmlpath",default = "./sample/imzml/", help = "imzml原始文件路径,默认./sample/imzml/")
  parser$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  parser$add_argument("-a","--asp",default = 1, type= "double",help = "长宽分辨率比")
  parser$add_argument("-m","--massrange",default = NULL,type= "double",help = "mz范围,默认NULL",
                      nargs = 2,dest = "mass.range")
  parser$add_argument("-re","--resolution",default = 100, type= "double",help = "分辨率,默认5")
  parser$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  
  parser$add_argument("-cf","--clusterfrom",default = "./sample/compress-raw/",help = "移动聚类结果路径,默认./sample/compress-raw/")
  parser$add_argument("-sp","--savepath",default = "./sample/compress-raw/", help = "imzml数据保存路径,默认./sample/compress-raw/")
  parser$add_argument("-bp","--areapath",default = "./sample/select/", help = "选区数据保存路径,默认./sample/select/")

  parser$add_argument("-ip","--infopath",default = "项目登记单.xlsx", help = "信息读取路径,默认项目登记单.xlsx/")
  args <- parser$parse_args()
  
  mulargs <- do.call(what = imzmlbgdata,args = args)
  
}
