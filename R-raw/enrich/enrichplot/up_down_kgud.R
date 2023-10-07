#!/opt/conda/bin/Rscript

#' 上下调对比图
#'
#' @export
up_down_kgud <- function(savepath = "./enrich/",  
			  incompare = "A_B", 
			  filt = "T", 
			  imagetype = c("pdf", "png") , 
			  height = 8, 
			  width = 14, 
			  dpi =300, 
			  fontfamily = "sans", 
                        ...){

    up_down(savepath = savepath, 
			  incompare = incompare, 
			  filt = filt, 
			  imagetype = imagetype , 
			  height = height, 
			  width = width, 
			  fontfamily = fontfamily, 
			  ...)
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
# 基本参数

  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 14, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 9, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")

  parser$add_argument("-ic","--incompare",default = "A_B", help = "比较组")
  parser$add_argument("-f","--filt",default = "T", help = "是否需要筛选")
  parser$add_argument("-s","--savepath",default = "./enrich/", help = "KEGG富集分析文件夹路径，默认./enrich/")

  args <- parser$parse_args()
  
  # 如果输入是压缩文件
  if (file.info(args$savepath)$isdir == FALSE) {
    
    if (dir.exists("enrich/KEGG")) {
      unlink("enrich/KEGG", recursive = TRUE)
    }
    dir.create("enrich/KEGG", showWarnings = FALSE, recursive = TRUE)
    system("chmod -R 777 enrich/KEGG")
	zip_filename <- tools::file_path_sans_ext(basename(args$savepath))
	args$incompare <- zip_filename
    unzip(args$savepath, exdir = "enrich/KEGG")
	#dir_list <- list.dirs("sample/", full.names = FALSE, recursive = FALSE)
	system(paste0("chmod -R 777 enrich/KEGG/",zip_filename))
	args$savepath <- "./enrich/"
  }
  
  result <- do.call(what = up_down_kgud,args = args)
}