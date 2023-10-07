#!/opt/conda/bin/Rscript

#' LC非靶数据获取
#'
#' @param data obj
#' @param ... 见`datafetch_untargetlcms`
#'
#' @export
datafetch_untargetlcms_obj <- function(data,
                                       ...){
  
  data$data$rawdata$data <- datafetch_untargetlcms(code.negID = data$info$datafrom$LC原始文件$negID,
                                                   code.posID = data$info$datafrom$LC原始文件$posID,
                                                   code.negM = data$info$datafrom$LC原始文件$negM,
                                                   code.posM = data$info$datafrom$LC原始文件$posM,
                                                   negname = data$info$sample$LCMS$neg$name,
                                                   posname = data$info$sample$LCMS$pos$name,
                                                   samplename =  data$info$sample$samplename,
                                                   weight =  data$info$sample$weight,
                                                   ...)
  
  return(data)
}

#' 获取QI定性定量文件
#' 
#' @param code.negID 负离子定性路径
#' @param code.posID 正离子定性路径 
#' @param code.negM 负离子M文件
#' @param code.posM 正离子M文件
#' @param negname 负离子M中样本名称 
#' @param posname 正离子M中样本名称 
#' @param samplename 样本分析名称
#' @param weight 样本称重
#' @param mode 处理模式，默认为"QI自带归一化模式"
#' @param deal mode为`raw`时,处理方式默认为"峰面积归一化"
#' @param neglab deal为`内标归一化`时,负离子内标名称
#' @param poslab deal为`内标归一化`时,正离子内标名称
#' @param data 定性数据
#' @param peptide 逻辑，是否保留肽段
#' @param meta 关注代谢物
#' @param score score筛选标准
#' @param fragscore 二级筛选标准
#' @param merge 逻辑值，是否仅保留未定性数据
#' @param ... 
#'
#' @export 
datafetch_untargetlcms <- function(code.negID,
                                   code.posID,
                                   code.negM,
                                   code.posM,
                                   negname = NULL,
                                   posname = NULL,
                                   samplename = NULL,
                                   weight = NULL,
                                   peptide = F,
                                   meta = NULL,
                                   score = 36,
                                   fragscore = 0,
                                   mode = "normal",
                                   deal = "峰面积归一化",
                                   neglab = NULL,
                                   poslab = NULL,
                                   merge = F,
                                   ...) {
  negposID <- dataID(code.negID = code.negID,
                     code.posID = code.posID,
                     peptide = peptide,
                     meta = meta,
                     score = score,
                     fragscore = fragscore,
                     ...)
  
  negposM <- dataM(code.negM = code.negM,
                   code.posM = code.posM,
                   negname = negname,
                   posname = posname,
                   samplename = samplename,
                   weight = weight,
                   mode = mode, deal = deal, neglab = neglab, poslab = poslab)
  
  data <- MIDmerge(dataM = negposM,
                   dataID = negposID,
                   merge = merge)
  
  return(data)
}
