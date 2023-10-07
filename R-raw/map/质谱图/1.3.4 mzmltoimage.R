#!/opt/conda/bin/Rscript

#' mzmap
#'
#' 将mzml数据生成质谱图
#'
#' @param path 数据路径
#' @param name 样本名称
#' @param deal 是否删除
#' @param type 输出格式
#' @param merge 是否叠加
#' @param ...
#'
#' @export
mzmap <- function(path,
                  name,
                  deal = NULL,
                  merge = F,
                  mergename = "Combine",
                  ...) {
  
  print("开始质谱转图")
  suppressMessages(library("MSnbase"))
  
  if (!is.null(deal)) {
    path <- path[!(!is.na(deal)&deal=="删除")]
    name <- name[!(!is.na(deal)&deal=="删除")]
  }
  
  if (any(is.na(path))) {
    name <- name[!is.na(path)]
    path <- path[!is.na(path)]
  }
  
  if(merge){
    
    map_mass_massmap(filename = path,
                     samplename = name,
                     mapname = mergename,
                     ...)
    
  }else{
    
    suppressWarnings(library(foreach))
    suppressWarnings(library(doParallel))
    registerDoParallel(cores=5)
    invisible(
    foreach(i=1:length(path)) %dopar% {
      # for (i in 1:length(path)) {
      
      map_mass_massmap(filename = path[i],
                       samplename = name[i],
                       mapname = name[i],
                       ...)
      
    })
  }
}

#' @export
mzmltoimage_obj <- function(filename = "内部分析单.xlsx",...){
  data <- readregistration(name = filename)
  mzmltoimage(data = data,...)
}

#' @export
mzmltoimage <- function(data,...) {
  UseMethod("mzmltoimage")
}

#' @export
mzmltoimage.default <- function(data,
                                ...) {
  print("~本项目类别不支持导图，跳过本步骤")
}


#' 非靶向LCMS生成质谱图
#'
#' @param data 项目信息
#' @param type 输出格式
#' @param ... 见[mzmap()]
#'
#' @export
mzmltoimage.UntargetLcms <- function(data,
                                     type = "BPC",
                                     ...) {
  
  mzmap(path = data$info$sample$LCMS$neg$mzml,
        name = paste0(data$info$sample$samplename,"-",type,"-LCMS-neg"),
        deal = data$info$sample$deal,
        type = type,
        ...)
  
  # mzmap(path = data$info$sample$LCMS$neg$mzml,
  #       name = paste0(data$info$sample$samplename,"-",type, "-LCMS-neg"),
  #       deal = data$info$sample$deal,
  #       type = type,
  #       merge = T,
  #       mergename = "Combine-LCMS-neg",
  #       ...)
  
  mzmap(path = data$info$sample$LCMS$pos$mzml,
        name = paste0(data$info$sample$samplename,"-",type, "-LCMS-pos"),
        deal = data$info$sample$deal,
        type = type,
        ...)
  
  # mzmap(path = data$info$sample$LCMS$pos$mzml,
  #       name = paste0(data$info$sample$samplename,"-",type, "-LCMS-pos"),
  #       deal = data$info$sample$deal,
  #       type = type,
  #       merge = T,
  #       mergename = "Combine-LCMS-pos",
  #       ...)
  
}

#' @export
mzmltoimage.UntargetLcmsEmdb <- mzmltoimage.UntargetLcms

#' @export
mzmltoimage.UntargetLcmsPmdb <- mzmltoimage.UntargetLcms

#' @export
mzmltoimage.UntargetLcmsTCM <- mzmltoimage.UntargetLcms

#' @export
mzmltoimage.UntargetLipid <- mzmltoimage.UntargetLcms

#' 非靶向GCMS生成质谱图
#'
#' @param data 项目信息
#' @param type 输出格式
#' @param ... 见[mzmap()]
#'
#' @export
mzmltoimage.UntargetGcms <- function(data,
                                     type = "BPC",
                                     ...) {
  
  mzmap(path = data$info$sample$GCMS$mzml,
        name = paste0(data$info$sample$samplename,"-",type, "-GCMS"),
        deal = data$info$sample$deal,
        type = type,
        ...)
  
}

#' @export
mzmltoimage.UntargetSpmeGcms <- mzmltoimage.UntargetGcms

#' @export
mzmltoimage.UntargetWaxGcms <- mzmltoimage.UntargetGcms


#' @export
mzmltoimage.UntargetBoth <- function(data,
                                     savepath = "./",
                                     ...){
  mzmltoimage.UntargetLcms(data = data,savepath = paste0(savepath,"/LCMS"),...)
  mzmltoimage.UntargetGcms(data = data,savepath = paste0(savepath,"/GCMS"),...)
}

#' @export
mzmltoimage.UntargetBothEmdb <- mzmltoimage.UntargetBoth

#' @export
mzmltoimage.UntargetBothPmdb <- mzmltoimage.UntargetBoth

#' @export
mzmltoimage.QtargetLipid <- function(data,
                                     type = "TIC",
                                     ...) {
  mzmltoimage.UntargetLcms(data = data,type = type,...)
}


#' @export
mzmltoimage.QtargetMic <- function(data,
                                   type = "TIC",
                                   ...) {
  mzmltoimage.UntargetBoth(data = data,type = type,...)
}
