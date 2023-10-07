#!/opt/conda/bin/Rscript

#' 提取空代背景后，剔除背景后压缩聚类
#'
#' @export
imzmlrmbgtocluster <- function(imzmlclusterpath = "./sample/compress/",
                               imzmlmovepath = "./sample/compress-raw/",
                               savepath = "./sample/compress/",
                               infopath = "项目登记单.xlsx",
                               bgimzmlpath = "./sample/compress-raw/",
                               bgintensity = 100,
                               bgfreq = 0.05,
                               multipleintensity = 2,
                               ...){
  
  
  if(file.exists(imzmlclusterpath)){
    if(!file.exists(imzmlmovepath)){
      file.rename(imzmlclusterpath,imzmlmovepath)
    }
  }
  
  if(file.exists(imzmlmovepath)){
    imzmlbgdata(clusterfrom = imzmlmovepath,
                savepath = bgimzmlpath,
                ...)
    AllCompressAndCluster(rebg = T,
                          bgimzmlpath = bgimzmlpath,
                          savepath = savepath,
                          bgintensity = bgintensity,
                          bgfreq =  bgfreq,
                          multipleintensity = multipleintensity,
                          ...)
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
  parser$add_argument("-re","--resolution",default = 100, type= "double",help = "分辨率,默认100")
  parser$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  
  parser$add_argument("-ic","--imzmlclusterpath",default = "./sample/compress/",help = "初始聚类结果路径,默认./sample/compress/")
  parser$add_argument("-im","--imzmlmovepath",default = "./sample/compress-raw/",help = "移动聚类结果路径,默认./sample/compress-raw/")
  
  parser$add_argument("-s","--savepath",default = "./sample/compress/", help = "imzml数据保存路径,默认./sample/compress/")
  parser$add_argument("-b","--bgimzmlpath",default = "./sample/compress-raw/", help = "背景imzml数据保存路径,默认./sample/compress-raw/")
  parser$add_argument("-bi","--bgintensity",default = 100,  type= "double",help = "背景离子识别强度")
  parser$add_argument("-bf","--bgfreq",default = 0.05,  type= "double",help = "背景离子表达范围")
  parser$add_argument("-mi","--multipleintensity",default = 2,  type= "double",help = "样本离子比背景离子表达倍数")
  
  parser$add_argument("-ip","--infopath",default = "项目登记单.xlsx", help = "信息读取路径,默认项目登记单.xlsx/")
  args <- parser$parse_args()
  
  mulargs <- do.call(what = imzmlrmbgtocluster,args = args)
  
}