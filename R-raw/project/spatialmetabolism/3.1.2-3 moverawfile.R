#!/opt/conda/bin/Rscript

#' 移动空代原始数据
#' 
#' @param imzmlpath 原始数据路径 
#' @param savepath 移动到路径
#'
#' @export
arrangerawfile <- function(imzmlpath = "./raw/",
                           savepath = "./sample/imzml/"){
  
  removeimzmlandmzml(path = imzmlpath)
  moverawfile(rawpath = imzmlpath,
              newpath = savepath)
  
}

#' @export
removeimzmlandmzml <- function(path = "raw"){
  
  if(dir.exists(paths = path)){
    removeibd <- list.files(path = path,pattern = "\\.ibd$",recursive = T,full.names = T)
    file.remove(removeibd)
    removeimzml <- list.files(path = path,pattern = "\\.imzML$",recursive = T,full.names = T)
    file.remove(removeimzml)
    removemzml <- list.files(path = path,pattern = "\\.mzML$",recursive = T,full.names = T)
    file.remove(removemzml)
    removecdf <- list.files(path = path,pattern = "\\.cdf$",recursive = T,full.names = T)
    file.remove(removecdf)
    removemsi <- list.files(path = path,pattern = "\\.msi$",recursive = T,full.names = T)
    file.remove(removemsi)
    removesld <- list.files(path = path,pattern = "\\.sld$",recursive = T,full.names = T)
    file.remove(removesld)
    removecache <- list.files(path = path,pattern = "\\.cache$",recursive = T,full.names = T)
    file.remove(removecache)
  }else{
    warning(paste0(path,"目录不存在"),immediate. = T)
  }

}

#' @export
moverawfile <- function(rawpath = "raw",
                        newpath = "./sample/imzml"){
  
  if(dir.exists(paths = rawpath)){
    system(paste0("mv -vf ",rawpath," ",newpath),
           ignore.stdout = F,
           ignore.stderr = F,
           wait = T)
  }else{
    warning(paste0(rawpath,"目录不存在"),immediate. = T)
  }
  
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-i","--imzmlpath",default = "./raw/", help = "imzml原始文件路径,默认./raw/")
  parser$add_argument("-s","--savepath",default = "./sample/imzml/", help = "imzml数据保存路径,默认./sample/imzml/")
  args <- parser$parse_args()
  
  writeinfo()
  
  mulargs <- do.call(what = arrangerawfile,args = args)
  
  writeinfo(endtime = T)
  
}
