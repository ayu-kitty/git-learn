#!/opt/conda/bin/Rscript

#' 对于拟靶向肠道定性
#'
#' @param obj 数据
#' @param methodpath 方法路径
#' @param outputpath 输出路径
#'
#' @export
QtargetQualitative_Mic_for_obj <- function(obj,
                                           outputpath = "./raw/搜库数据",
                                           methodpath = "database/qualitative/QtargetMic"){
  data <- obj
  mzmlData <- QtargetQualitative(filename = data$info$sample$LCMS$neg$mzml,
                                 methodpath = selectfile(path = methodpath,file = "qualitative-neg.txt"),
                                 outputpath = outputpath,
                                 mode = "neg")
  data$info$datafrom$LC原始文件$negID <- mzmlData[["ouputpath"]]
  mzmlData <- QtargetQualitative(filename = data$info$sample$LCMS$pos$mzml,
                                 methodpath = selectfile(path = methodpath,file = "qualitative-pos.txt"),
                                 outputpath = outputpath,
                                 mode = "pos")
  data$info$datafrom$LC原始文件$posID <- mzmlData[["ouputpath"]]
  
  return(data)
}
