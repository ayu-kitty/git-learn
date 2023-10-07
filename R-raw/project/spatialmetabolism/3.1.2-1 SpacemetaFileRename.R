#!/opt/conda/bin/Rscript

#' 空间代谢组原始数据名检查
#'
#' @param datatype 数据类型
#' @param path 检查路径
#' @param ...
#'
#' @export
SpacemetaFileRename <- function(datatype = "\\.raw$",
                                path = "raw",
                                ...) {
  # datatype = "\\.raw$"
  # path = "raw"
  imzMLdir <- list.dirs(path = path, full.names = T)[-1]
  imzMLdir2 <- list.dirs(path = path, full.names = F)[-1]
  
  if(length(imzMLdir) == 0){
    return()
  }
  
  for (i in seq_len(length(imzMLdir))) {
    imzMLfile <- list.files(path = imzMLdir[i], pattern = datatype, full.names = T)
    imzMLfile2 <- list.files(path = imzMLdir[i], pattern = datatype, full.names = F)
    
    # 判断目录下是否有imzML文件
    if (length(imzMLfile) == 0) {
      
    } else {
      if (all(nchar(imzMLfile2)[1] == nchar(imzMLfile2))) {
        print(paste0(imzMLdir[i], "目录无问题不进行修改"))
      } else {
        print(paste0(imzMLdir[i], "目录有问题进行修改"))
        
        imzMLfile3 <- imzMLfile2
        # 原始数据名称修改,补0,统一位数
        for (j in 2:max(nchar(imzMLfile3))) {
          if (all(nchar(imzMLfile3)[1] == nchar(imzMLfile3))) {
            break
          } else if (all(substring(imzMLfile2, j, j)[1] == substring(imzMLfile2, j, j))) {
            next
          } else {
            imzMLfile3[nchar(imzMLfile3) != max(nchar(imzMLfile3))] <- gsub(
              pattern = substring(imzMLfile2[1], 1, j - 1),
              replacement = paste0(substring(imzMLfile2[1], 1, j - 1), "0"),
              x = imzMLfile3[nchar(imzMLfile3) != max(nchar(imzMLfile3))]
            )
          }
        }
        
        file.rename(paste0(imzMLdir[i], "/", imzMLfile2), paste0(imzMLdir[i], "/", imzMLfile3))
      }
    }
  }
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-p","--path",default = "./raw/", help = "原始文件路径,默认./raw/")
  args <- parser$parse_args()
  
  writeinfo()
  
  mulargs <- do.call(what = SpacemetaFileRename,args = args)
  
  writeinfo(endtime = T)

}
