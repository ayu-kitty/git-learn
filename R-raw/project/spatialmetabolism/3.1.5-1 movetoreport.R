#!/opt/conda/bin/Rscript

#' movetoreport
#'
#' 生成空代项目文件
#' samplepath sample目录路径
#' type 空代产品线名称(空间代谢组;中药空间代谢组)
#' @export

movetoreport <- function(asp = 1,
                         samplepath = "../sample/") {
  ptm <- Sys.time()
  
  report <- getreportname()
  
  setwddir(filename = report, force = T)
  
  copydir(from = "../HE/",
          to = "./1.成像图/HE染色/")
  
  Mulgetallrdsdata(rdspath = paste0(samplepath,"cluster/"),
                   savepath = "1.成像图/SSCC/",
                   qualitativepath = paste0(samplepath,"qualitative/Qualitative.xlsx"),
                   asp = asp)
  
  if(length(list.files(path = paste0(samplepath,"/final/"),pattern = "^qc_data")) > 0){
    imzmlpath <- paste0(samplepath,"/adjustdata/")
  }else{
    imzmlpath <- paste0(samplepath,"/final/")
  }
  
  Mulplotmap(imzmlpath = imzmlpath,
             savepath = "1.成像图/质谱图")
  Mulimagemap(imzmlpath = imzmlpath,
              savepath = "1.成像图/成像图/",
              mapmz = T,
              vague = F,
              asp = asp)
  # Mulimagemap(imzmlpath = paste0(samplepath,"final/"),
              # savepath = "1.成像图/成像图-平滑化/",
              # mapmz = T,
              # vague = T,
              # asp = asp)
  
  copydir(from = paste0(samplepath,"qualitative/Qualitative.xlsx"),
          to = "./2.定性结果/")
  
  getareadatainreport(qualitativefile = paste0(samplepath,"qualitative/Qualitative.xlsx"),
                      areadatapath = paste0(samplepath,"area/data/"),
                      savepath = "./3.选区数据/数据矩阵")
  
  copydir(from = paste0(samplepath,"area/map/"),
          to = "./3.选区数据/成像图")
 
 # 根据是否有比较分析修改文件名
  if (length(list.files(paste0(samplepath,"analysis")))> 0) {
	  if(dir.exists(paste0(samplepath,"analysis/raw"))){
		copydir(from = paste0(samplepath,"analysis/result/"),
				to = "./4.比较分析/")
	  }else{
		copydir(from = paste0(samplepath,"analysis/"),
				to = "./4.比较分析/",
				rmfile = "raw.RData")
	  }
  copydir(from = databasepath(path = "script/spatialmetabolism/区域选择/"),
          to = "./5.原始数据/")
  copydir(from = paste0(samplepath,"pre/"),
          to = "./5.原始数据/原始数据/",
          rmfile = "^qc_data")
  copydir(from = "../样品登记单.xlsx",
          to = "./6.样品登记单/")
} else{
  copydir(from = databasepath(path = "script/spatialmetabolism/区域选择/"),
          to = "./4.原始数据/")
  copydir(from = paste0(samplepath,"pre/"),
          to = "./4.原始数据/原始数据/",
          rmfile = "^qc_data")
  copydir(from = "../样品登记单.xlsx",
          to = "./5.样品登记单/")
  
}
	# 复制成像图，SSCC和质谱图到176盘进行检查
	pattern <- "DZLM(\\d+)"
	AnalysisId <- regmatches(report, gregexpr(pattern, report))
	path_176 <- paste0("/data/nas/176/空代-空代实验组-jpp/项目数据/2023/",AnalysisId)
	if (!dir.exists(path_176)) {
	  dir.create(path_176)}
	copydir(from = "./1.成像图/成像图/",
			to = paste0(path_176,"/成像图") )
	copydir(from = "./1.成像图/质谱图/",
			to = paste0(path_176,"/质谱图") )
	copydir(from = "./1.成像图/SSCC/",
			to = paste0(path_176,"/SSCC") )
	copydir(from = "./2.定性结果/",
			to = paste0(path_176) )
			
  
  setwd("../")
  
  ptm1 <- Sys.time() - ptm
  print(paste0("文件转移总共花费时间为", ptm1 / 60, "分钟"))
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  parser$add_argument("-a","--asp",default = 1, type= "double",help = "长宽分辨率比")
  parser$add_argument("-t","--type",default = "空间代谢组",help = "空代产品线类型:空间代谢组 or 中药空间代谢组 or 空间代谢组-waters")
  parser$add_argument("-s","--samplefrom",default = "./sample/final/",help = "样本数据路径,默认./sample/final/")
  parser$add_argument("-sp","--samplepath",default = "../sample/",help = "sample目录路径,默认 ../sample/")
  
  args <- parser$parse_args()
  
  writeinfo()
  
  movetoreport(asp = args$asp,samplepath = args$samplepath)
  makespacereport(type = args$type)
  
  writeinfo(endtime = T)
  
}
