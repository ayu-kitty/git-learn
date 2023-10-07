#!/opt/conda/bin/Rscript

#' 保存imzml文件
#'
#' @param data 数据
#' @param filename 文件名
#'
#' @export
saveimzml <- function(data,
                      filename,
                      overwrite = T,
                      ...) {
  suppressMessages(library("Cardinal"))
  if(is.null(data)) {return()}
  
  print(paste0("在", filename, "中保存"))
  
  data <- data[, order(coord(data)$y, coord(data)$x)]
  
  if(!file.exists(filename) | overwrite){
    filename <- imzml_unlink(savepath = dirname(filename),
                             samplename = basename(gsub(pattern = "\\.imzML",
                                                        replacement = "",
                                                        x = basename(filename))))
    rdsfilename <- gsub(pattern = "\\.imzML",replacement = ".ibr",x = filename)
    ibdfilename <- gsub(pattern = "\\.imzML",replacement = ".ibd",x = filename)
    
    writeMSIData(object = data,
                 file = filename,
                 outformat = "imzML",
                 mz.type = "64-bit float",
                 ...)
    
    ibrdata <- list(pData = pData(data))
    saverds(data = ibrdata,filename = rdsfilename)
    
    Sys.chmod(paths = filename,mode = "0777",use_umask = F)
    Sys.chmod(paths = ibdfilename,mode = "0777",use_umask = F)
    Sys.chmod(paths = rdsfilename,mode = "0777",use_umask = F)
  }

  # print("保存完毕")
}

#' @export
imzml_unlink <- function(savepath,
                         samplename){
  
  # 判断slide是否存在，不存在进行创建
  savepath <- createdir(filename = savepath)
  
  filename <- paste0(savepath, "/",samplename, ".imzML")
  filename_ibd <- paste0(savepath, "/",samplename, ".ibd")
  filename_ibr <- paste0(savepath, "/",samplename, ".ibr")
  
  # 判断文件是否存在，存在进行删除
  if (file.exists(filename) | file.exists(filename_ibd) | file.exists(filename_ibr)) {
    unlink(filename)
    unlink(filename_ibd)
    unlink(filename_ibr)
  }
  
  return(filename)
}
