
#' datafetch
#'
#' 数据获取
#'
#' @param data obj
#'
#' @export
datafetch <- function(data,
                      ...) {
  UseMethod("datafetch")
}

#' @export
datafetch.default <- function(data,
                              ...) {
  stop("项目类别未选择或填错")
}


#' LC非靶数据获取
#'
#' @export
datafetch.UntargetLcms <- function(data,
                                   rttimemin = 0.5,
                                   rttimemax = 15.5,
                                   ...) {
  
  data <- datafetch_untargetlcms_obj_2(data = data,...)
  
  data <- rentiontimefilter(data = data,rttimemin=rttimemin,rttimemax=rttimemax)
  
  return(data)
}

#' LC-EMDB非靶数据获取
#'
#' @export
datafetch.UntargetLcmsEmdb <- function(data,
                                       ...) {
  
  data <- datafetch.UntargetLcms(data = data,...)
  
  return(data)
}

#' LC-PMDB非靶数据获取
#'
#' @export
datafetch.UntargetLcmsPmdb <- function(data,
                                       ...) {
  
  data <- datafetch.UntargetLcms(data = data,...)
  
  return(data)
}

#' 中药非靶数据获取
#'
#' @export
datafetch.UntargetLcmsTCM <- function(data,
                                      Tcm = T,
                                      ...) {
  
  data <- datafetch.UntargetLcms(data = data,Tcm = Tcm,...)
  
  return(data)
}

#' GC非靶数据获取
#'
#' @param data obj
#' @param mode area峰面积,fame内标分段归一化,single单一内标归一化
#' @param famersd 0.1 内标rsd筛选标准0.1
#' @param projectclass spme(顶空) or gc(非靶gc) or wax(蜡质)
#' @param RI T of F 是否进行RI计算
#' @param score 筛选打分默认70
#' @param versions 3:lug3.0;4:lug4.0(2022.10.24更新)
#' @param ... 见[gcpredeal()]
#'
#' @export
datafetch.UntargetGcms <- function(data,
                                   mode = "fame",
                                   RI = T,
                                   score = 70,
                                   projectclass = "gc",
                                   famersd = 0.1,
                                   ...) {
  
  data <- datafetch_untargetgcms_obj_2(data = data,
                                       mode = mode,RI = RI,
                                       samplename = data$info$sample$samplename,
                                       grouping = data$info$sample$rawname,
                                       score = score,projectclass = projectclass,famersd = famersd,
                                       ...)
  
  return(data)
}


#' GC非靶顶空数据获取
#'
#' @param data obj
#' @param mode area峰面积,fame内标分段归一化,single单一内标归一化
#' @param famersd 0.1 内标rsd筛选标准0.1
#' @param projectclass spme(顶空) or gc(非靶gc) or wax(蜡质)
#' @param RI T of F 是否进行RI计算
#' @param score 筛选打分默认70
#' @param versions 3:lug3.0;4:lug4.0(2022.10.24更新)
#' @param ... 见[gcpredeal()]
#'
#' @export
datafetch.UntargetSpmeGcms <- function(data,
                                       mode = "area",
                                       RI = F,
                                       score = 70,
                                       projectclass = "spme",
                                       famersd = 0.1,
                                       versions = 4,
                                       ...) {
  
  data <- datafetch_untargetgcms_obj_2(data = data,
                                       mode = mode,RI = RI,
                                       score = score,projectclass = projectclass,famersd = famersd,
                                       ...)
  
  return(data)
}

#' GC非靶蜡质数据获取
#'
#' @param data obj
#' @param mode area峰面积,fame内标分段归一化,single单一内标归一化
#' @param famersd 0.1 内标rsd筛选标准0.1
#' @param projectclass spme(顶空) or gc(非靶gc) or wax(蜡质)
#' @param RI T of F 是否进行RI计算
#' @param score 筛选打分默认70
#' @param versions 3:lug3.0;4:lug4.0(2022.10.24更新)
#' @param ... 见[gcpredeal()]
#'
#' @export
datafetch.UntargetWaxGcms <- function(data,
                                      mode = "area",
                                      RI = F,
                                      score = 70,
                                      projectclass = "wax",
                                      famersd = 0.1,
                                      versions = 4,
                                      ...) {
  
  data <- datafetch_untargetgcms_obj_2(data = data,
                                       mode = mode,RI = RI,
                                       score = score,projectclass = projectclass,famersd = famersd,
                                       ...)
  
  return(data)
}


#' 外来数据获取 
#' @export
datafetch.UntargetOut <- function(data,
                                  filename = data$info$datafrom$其他原始文件,
                                  ...) {
  
  data$data$rawdata$data <- readdata(filename = filename, ...)
  
  return(data)
}

#' @export
datafetch.OutData <- function(data,
                              sheet = "数据矩阵",
                              ...) {
  
  data <- datafetch.UntargetOut(data,
                                sheet = sheet,
                                ...)
  
  return(data)
}

#' 拟靶向数据获取
#'
#' @param data obj
#' @param methodfrom 定性方法文件
#' @param methodpath 定性方法路径
#' @param ... 见[datafetch.QtargetLipid()]
#'
#' @export
datafetch.QtargetLipid <- function(data,
                                   ...) {
  
  if (!is.na(data$info$datafrom$LC原始文件$negID) |
      !is.na(data$info$datafrom$LC原始文件$posID)) {
    
    data$data$rawdata$data <- QTargetedLipidDataM(data, ...)
    
    return(data)
  } else {
    data <- QtargetQualitative_Lipid_for_obj(data)
    
    data <- datafetch.QtargetLipid(data,
                                   ...)
    return(data)
  }
}

#' 肠道拟靶向数据获取
#'
#' @param data obj
#' @param ...
#'
#' @export
datafetch.QtargetMic <- function(data,
                                 ...) {
  
  if (!is.na(data$info$datafrom$LC原始文件$negID) |
      !is.na(data$info$datafrom$LC原始文件$posID)) {
    
    if(is.na(data$info$datafrom$GC原始文件)){
      gcfile <- list.files(path = "./",pattern = "^GCMS")
      if(length(gcfile) > 0){
        data$info$datafrom$GC原始文件 <- gcfile
      }else{
        stop("未找到GCMS定性文件")
      }
    }
    
    data <- micfetch(data, ...)
    data$info$basic$索引$`Ion mode` <- NA
    
    return(data)
  } else {
    
    data <- QtargetQualitative_Mic_for_obj(data)
    
    data <- datafetch.QtargetMic(data,
                                 ...)
    return(data)
  }
  
  return(data)
}


#' @export
datafetch.Untarget4DLipid <- function(data,
                                      sheet = "数据矩阵",
                                      ...) {
  data <- datafetch.UntargetOut(data,
                                sheet = sheet,
                                ...)
  
  return(data)
}

#' 非靶向代谢-脂质数据获取
#'
#' @param data obj
#' @param ...
#'
#' @export
datafetch.UntargetLipid <- function(data,
                                    ...) {
  
  data$data$rawdata$data <- lipid_predeal(filename_neg=data$info$datafrom$LC原始文件$negID,
                                          filename_pos=data$info$datafrom$LC原始文件$posID,
                                          samplename=data.frame(rawname = data$info$sample$rawname,
                                                                samplename = data$info$sample$samplename,
                                                                negname = data$info$sample$LCMS$neg$name,
                                                                posname = data$info$sample$LCMS$pos$name,
                                                                stringsAsFactors = F))
  
  return(data)
}

#' 靶向代谢-胆汁酸数据获取
#'
#' @param inputpath 分析数据输入路径
#' @param tlcms TLCMS-胆汁酸.xlsx
#' @param qcinfo 标品信息表.xlsx
#' @param rawData 下机数据
#'
#' @export
datafetch.Target <- function(data,
                             ...){
  
  rawdata_file <- list.files(pattern = "*-results-*")
  tlcms_file <- list.files(pattern = "TLCMS-.*\\.xlsx")
  qcinfo <- list.files(pattern = "*-标品信息表.xlsx")
  data$data$rawdata$data <- datafech_targettlcms (
                              tlcms = tlcms_file,
                              qcinfo = qcinfo,
                              rawData = rawdata_file)
  return(data)
}
