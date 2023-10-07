#!/opt/conda/bin/Rscript


#' 获取登记单
#'
#' @export
Getspaceregistration <- function(filepath = databasepath(database = "database/script/spatialmetabolism/"),
                                 overwrite = F) {
  print("项目登记单.xlsx复制中")
  file.copy(
    from = paste0(filepath,"/项目登记单.xlsx"),
    to = "./", overwrite = overwrite
  )

}


#' 获取背景扣除脚本
#'
#' @export
GetRejectBGscript <- function(filepath = databasepath(database = "database/script/spatialmetabolism/"),
                              overwrite = F) {
  print("RejectBG-Manual.R复制中")
  file.copy(
    from = paste0(filepath,"/RejectBG-Manual.R"),
    to = "./", overwrite = overwrite
  )
  
  print("SelectArea.R复制中")
  file.copy(
    from = paste0(filepath,"/SelectArea.R"),
    to = "./", overwrite = overwrite
  )

}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filepath",default = "/luming/script/spatialmetabolism/", help = "基础文件复制路径")
  parser$add_argument("-o","--overwrite",default = F, help = "是否覆盖原文件",action='store_true')
  args <- parser$parse_args()
  
  writeinfo()
  
  createdir(filename = "./sample/qualitative",linkdir = T)
  createdir(filename = "./sample/select",linkdir = T)
  
  result <- Getspaceregistration(filepath = args$filepath,
                                 overwrite = args$overwrite)
  result <- GetRejectBGscript(filepath = args$filepath,
                              overwrite = args$overwrite)
  
  writeinfo(endtime = T)
  
}
