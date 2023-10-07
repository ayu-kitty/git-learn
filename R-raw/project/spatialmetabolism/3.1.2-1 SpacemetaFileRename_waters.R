#!/opt/conda/bin/Rscript

#' 空间代谢组原始数据名检查
#'
#' @param datatype 数据类型
#' @param path 检查路径
#' @param ...
#'
#' @export
SpacemetaFileRename_waters <- function(datatype = "\\.imzml$",
                                       path = "raw",
                                       ...) {
  imzMLfile <- list.files(path = path,pattern = datatype,full.names = T,recursive = T)
  
  if (length(imzMLfile) > 0) {
    
    # 修改原始数据名称从.imzml为.imzML后缀
    changeimzMLfile <- gsub(pattern = datatype,replacement = ".imzML",x = imzMLfile)
    
    for ( i in 1:length(imzMLfile)) {
      
      file.rename(from = imzMLfile[i],to = changeimzMLfile[i])
      
    }
    
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-p","--path",default = "./raw/", help = "原始文件路径,默认./raw/")
  args <- parser$parse_args()
  
  writeinfo()
  
  mulargs <- do.call(what = SpacemetaFileRename_waters,args = args)
  
  writeinfo(endtime = T)
  
}
