#!/opt/conda/bin/Rscript

#' 降维聚类的运算和绘图
#'
#' @export
imzmltosscc_mulgetdata <- function(imzmlpath = "./sample/final/",
                        radius = 1,
                        clusters = 8,
                        sparse = 3,
                        mulk = 15,
						muls = 3,
						mulr =  1,
						allanalyst = "T",
                        ...){

    MulanaCluster(imzmlpath = imzmlpath,
                          r = radius,
                          k = clusters,
                          s = sparse,
                          mulk = mulk,
                          muls = muls,
                          mulr = mulr,
                          allanalyst = allanalyst,
						  ...)
    Mulgetallrdsdata() 

}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-i","--imzmlpath",default = "./sample/final/", help = "输入数据目录或者压缩文件")
  parser$add_argument("-k","--clusters",default = 8,type= "integer", help = "聚类数")
  parser$add_argument("-s","--sparse",default = 3,type= "integer", help = "稀疏参数")
  parser$add_argument("-r","--radius",default = 1,type= "integer", help = "距离半径")
  parser$add_argument("-mk","--mulk",default = 15,type= "integer", help = "多样本聚类数")
  parser$add_argument("-ms","--muls",default = 3,type= "integer", help = "多样本稀疏参数")
  parser$add_argument("-mr","--mulr",default = 1,type= "integer", help = "多样本距离半径")
  parser$add_argument("-al","--allanalyst",default = "T",help = "是否对all进行sscc")

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
	args$imzmlpath <- paste0("sample/",dir_list,"/")
	system(paste0("chmod -R 777 ",args$imzmlpath))
	dir.create("sample/map", showWarnings = FALSE, recursive = TRUE)
	dir.create("sample/cluster", showWarnings = FALSE, recursive = TRUE)
	system("chmod -R 777 sample/{map,cluster}")
  }
  #createdir(filename = "./sample/cluster",linkdir = T)
  args$allanalyst <- as.logical(args$allanalyst)
  
  result <- do.call(what = imzmltosscc_mulgetdata,args = args)
}