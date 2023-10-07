#!/opt/conda/bin/Rscript

#' movedata
#'
#' 上传obs
#'
#' @param data 项目信息
#' @param obsyear 存储位置
#' @param path 上传路径
#' @param ...
#'
#' @export

moverawdata <- function(rawpath,
                        deal,
                        to) {
  rawpath2 <- rawpath[!(deal %in% "删除")]
  rawpath2 <- rawpath2[!is.na(rawpath2)]
  if (length(rawpath2) == 0) {
    warning("无原始数据转移")
  } else {
    meta::dirfile(to)
    if (grepl(pattern = ".wiff$", x = rawpath2[1])) {
      rawpath2 <- c(rawpath2, paste0(rawpath2, ".scan"))
    }
    file.copy(from = rawpath2, to = to, recursive = T, overwrite = T)
  }
}


#' @export
movedata <- function(data) {
  UseMethod("movedata")
}

#' @export
movedata.default <- function(data) {
  warning(paste0(class(data), "项目类别不进行原始数据转移"))
}

#' @export
movedata.UntargetLcms <- function(data,
                                  path = "./raw/质谱数据/发送数据/") {
  print("LCMS负离子数据转移")
  moverawdata(
    rawpath = data$info$sample$LCMS$neg$raw,
    deal = data$info$sample$deal,
    to = paste0(path, "LCMS/neg")
  )
  print("LCMS正离子数据转移")
  moverawdata(
    rawpath = data$info$sample$LCMS$pos$raw,
    deal = data$info$sample$deal,
    to = paste0(path, "LCMS/pos")
  )

  # system(paste0("cd ", path, ";", databasepath(path = "script/md5/SingleFileMD5.sh")))
}

#' @export
movedata.UntargetLcmsEmdb <- function(data,
                                      ...) {
  movedata.UntargetLcms(
    data = data,
    ...
  )
}

#' @export
movedata.UntargetLcmsPmdb <- function(data,
                                      ...) {
  movedata.UntargetLcms(
    data = data,
    ...
  )
}
#' @export
movedata.UntargetLcmsTCM <- function(data,
                                      ...) {
  movedata.UntargetLcms(
    data = data,
    ...
  )
}
#' @export
movedata.UntargetLipid <- function(data,
                                  path = "./raw/质谱数据/发送数据/") {
  print("脂质负离子数据转移")
  moverawdata(
    rawpath = data$info$sample$LCMS$neg$raw,
    deal = data$info$sample$deal,
    to = paste0(path, "脂质/neg")
  )
  print("脂质正离子数据转移")
  moverawdata(
    rawpath = data$info$sample$LCMS$pos$raw,
    deal = data$info$sample$deal,
    to = paste0(path, "脂质/pos")
  )

  # system(paste0("cd ", path, ";", databasepath(path = "script/md5/SingleFileMD5.sh")))
}

#' @export
movedata.QtargetLipid <- function(data,
                                   path = "./raw/质谱数据/发送数据/") {
  print("拟靶向脂质负离子数据转移")
  moverawdata(
    rawpath = data$info$sample$LCMS$neg$raw,
    deal = data$info$sample$deal,
    to = paste0(path, "拟靶向脂质/neg")
  )
  print("拟靶向脂质正离子数据转移")
  moverawdata(
    rawpath = data$info$sample$LCMS$pos$raw,
    deal = data$info$sample$deal,
    to = paste0(path, "拟靶向脂质/pos")
  )

  # system(paste0("cd ", path, ";", databasepath(path = "script/md5/SingleFileMD5.sh")))
}



#' @export
movedata.UntargetGcms <- function(data,
                                  path = "./raw/质谱数据/发送数据/") {
  print("GCMS数据转移")
  moverawdata(
    rawpath = data$info$sample$GCMS$raw,
    deal = data$info$sample$deal,
    to = paste0(path, "GCMS")
  )
  # system(paste0("cd ", path, ";", databasepath(path = "script/md5/SingleFileMD5.sh")))
}

#' @export
movedata.UntargetSpmeGcms <- function(data,
                                      path = "./raw/质谱数据/发送数据/") {
  print("GCMS顶空数据转移")
  moverawdata(
    rawpath = data$info$sample$GCMS$raw,
    deal = data$info$sample$deal,
    to = paste0(path, "GCMS-Spme")
  )
  # system(paste0("cd ", path, ";", databasepath(path = "script/md5/SingleFileMD5.sh")))
}


#' @export
movedata.UntargetWaxGcms <- function(data,
                                     path = "./raw/质谱数据/发送数据/") {
  print("GCMS蜡质数据转移")
  moverawdata(
    rawpath = data$info$sample$GCMS$raw,
    deal = data$info$sample$deal,
    to = paste0(path, "GCMS-Wax")
  )
  # system(paste0("cd ", path, ";", databasepath(path = "script/md5/SingleFileMD5.sh")))
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-d","--data",default = "../data", help = "data.RData")
  
  args <- parser$parse_args()
  
  result <- do.call(what = movedata,args = args) 
  
}

