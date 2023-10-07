#!/opt/conda/bin/Rscript

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-r","--ratio",default = 2, type= "double",help = "仪器每秒扫描数量,默认2",required = T)
  parser$add_argument("-m","--massrange",default = NULL,type= "double",help = "mz范围",
                      nargs = 2,dest = "mass.range")
  parser$add_argument("-re","--resolution",default = 5, type= "double",help = "分辨率,默认5")
  parser$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  parser$add_argument("-a","--asp",default = 1, type= "double",help = "长宽比")
  parser$add_argument("-w","--waters",default = F,help = "是否使用waters仪器",action ='store_true')
  parser$add_argument("-tr","--tolerance",default = 0, type= "double",help = "waters峰校正参数,默认20")
  parser$add_argument("-c","--code",help = "执行编号", required = T)
  parser$add_argument("-id","--AnalysisId",default = "none",help = "项目编号")
  args <- parser$parse_args()
  
  writeinfo()
  
  # 创建基础脚本及登记单
  Getspaceregistration()
  GetRejectBGscript()
  # 自动填写登记单
  system(paste0("source /etc/profile;source ~/.bashrc;",  # 执行系统命令
                packagepath(path = "python/lmbio/spatialmetabolism/GetRegistrationInfo.py"),  # 使用 packagepath 获取指定路径
                " -c ",args$code,  # 设置参数 
                " -id ",args$AnalysisId
  ))
  
  createdir(filename = "./sample/imzml",linkdir = T)
  createdir(filename = "./sample/qualitative",linkdir = T)
  createdir(filename = "./sample/select",linkdir = T)
  
  if(!args$waters){
    # 非waters 一起数据转格式
    args$waters <- NULL
    args$tolerance <- NULL
    # 修改文件名,raw文件名称补齐2-3位,例如:从1改为01
    SpacemetaFileRename()
    # 转格式
    system(packagepath("command/spatialmetabolism_1.2_convertrawtoimzml.sh"))
    
    # imzml文件合并,按行合并为一个raw
    mulargs <- do.call(what = AllMulImzMLToOne,args = args)
    # 移动原始文件 从raw转移到/sample/imzml
    arrangerawfile()
  }else{
    args$waters <- NULL
    args$ratio <- NULL
	
	# 创建waterspeak文件夹，并复制表格
	if (file.exists("./excel") && file.info("./excel")$isdir) {
	  if (length(list.files("excel")) > 0) {
		createdir(filename = "./sample/waterspeak")
		excel_directory <- "./excel/"
		waterspeak_path <- "./sample/waterspeak/"
		excel_files <- list.files(excel_directory, pattern = "\\.csv$", ignore.case = TRUE)
		for (file in excel_files) {
		  if (grepl("neg|pos", file, ignore.case = TRUE)) {
			source_path <- file.path(excel_directory, file)
			file.copy(source_path, waterspeak_path)
			  }
			}
		  }
		 }
		  
    # 运行 waters 文件名修正
    SpacemetaFileRename_waters()
    # water数据处理,包括绘制内标mz556,mz554的质谱图存入imzml文件夹
    mulargs <- do.call(what = imzmlwatersdatadeal,args = args)
  }
  
  
  # 压缩及聚类
  args$ratio <- NULL
  args$resolution <- 100
  args$units <- "ppm"
  # 压缩聚类用于选区,保存于./sample/compress/
  mulargs <- do.call(what = AllCompressAndCluster,args = args)
  
  writeinfo(endtime = T)
  
}
