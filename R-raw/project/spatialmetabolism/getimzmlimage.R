#!/opt/conda/bin/Rscript

#' 批量绘制空代质谱成像图成像图
#'
#' @export
get_imzmlimage <- function(imzmlpath = "./sample/final/",
                           imagetype = c("jpg", "pdf"),
						    width = 8,
                            height = 7,
							color = c("blue2", "cyan2", "yellow",
                                 "brown1", "firebrick3"),
							resolution = 5,
							xlab = "",
                            ylab = "",
							...
    ){
   
   #批量绘制累加质谱成像图
   getsumimzmlimage(imzmlpath= imzmlpath,
                            imagetype = imagetype,
						    width = width,
                            height = height,
							color = color,
							resolution = resolution,
							xlab = xlab,
                            ylab = ylab,
							...
   
   )

   #批量绘制表达物质数的质谱成像图
   getnumberimzmlimage(imzmlpath= imzmlpath,
                            imagetype = imagetype,
						    width = width,
                            height = height,
							color = color,
							resolution = resolution,
							xlab = xlab,
                            ylab = ylab,
							...
   )
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  suppressMessages(library("meta"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-i","--imzmlpath",default = "./sample/final/", help = "输入数据目录或者压缩文件")
  parser$add_argument("-im","--imagetype",default = c("jpg", "pdf"), help = "图片格式")
  parser$add_argument("-w","--width",default =8,type= "double", help = "图片宽度")
  parser$add_argument("-g","--height",default = 7, type= "double", help = "图片高度")
  parser$add_argument("-c","--color",default =c("blue2", "cyan2", "yellow","brown1", "firebrick3") , help = "成像颜色")
  parser$add_argument("-r","--resolution",default = 5,type= "double", help = "图像分辨率")
  parser$add_argument("-x","--xlab",default = "", help = "x轴标签")
  parser$add_argument("-y","--ylab",default = "", help = "y轴标签")

  args <- parser$parse_args()
  # 如果输入是压缩文件
  if (file.info(args$imzmlpath)$isdir == FALSE) {
    if (dir.exists("sample")) {
      unlink("sample", recursive = TRUE)
    }
    dir.create("sample/", showWarnings = FALSE, recursive = TRUE)
	system("chmod -R 777 sample")
    unzip(args$imzmlpath, exdir = "sample/")
	dir_list <- list.dirs("sample/", full.names = FALSE, recursive = FALSE)
	args$imzmlpath <- paste0("./sample/",dir_list,"/")
	system(paste0("chmod -R 777 ",args$imzmlpath))
	dir.create("sample/map/Intensity", showWarnings = FALSE, recursive = TRUE)
	system("chmod -R 777 sample/map/*")
  }
  result <- do.call(what = get_imzmlimage,args = args)
}
