#!/opt/conda/bin/Rscript

#' 根据meta包中的obj保存多元统计的yaml文件
#' 
#' @param obj 数据
#' @param datafile 数据保存路径
#' @param classfile 分组保存路径
#' @param mulstatisticsyamlfile 多元统计保存路径
#' @param mulstatisticsyamlfile 多元统计yaml保存路径
#'
#' @export
writemulstatisticsyaml <- function(obj,
                                   datafile = "oecloud/rawdata/datafile.txt",
                                   classfile = "oecloud/rawdata/classfile.yaml",
                                   mulstatistics_rds_savepath = "oecloud/mulstatisticsanalyst",
                                   mulstatisticsyamlfile = "oecloud/rawdata/mulstatistics.yaml",
                                   all = T){
  data <- obj
  compare <- list(compare = list(),
                  path = list(datafile = datafile,
                              classfile = classfile,
                              mulstatistics_rds_savepath = mulstatistics_rds_savepath))
  for (i in 1:length(data$info$dif$compare)) {
    for (j in 1:length(data$info$dif$compare[[i]])) {
      comparename <- gsub(pattern = "/",replacement = "-vs-",x = data$info$dif$compare[[i]][[j]])
      group <- unlist(strsplit(x = data$info$dif$compare[[i]][[j]],split = "/"))
      compare$compare[[comparename]] <- list(rawname = data$info$dif$compare[[i]][[j]],
                                             group = group,
                                             "PCA" = list(scaleC = data$info$dif$PCAscaleC[[i]][j],
                                                          predI = -1,
                                                          orthoI = 0,
                                                          permI = 0,
                                                          log10L = data$info$dif$log10L[[i]][j],
                                                          amount = 3),
                                             "PLS-DA" = list(scaleC = data$info$dif$PLSscaleC[[i]][j],
                                                             predI = -1,
                                                             orthoI = 0,
                                                             permI = ifelse(length(group) > 2,200,0),
                                                             log10L = data$info$dif$log10L[[i]][j],
                                                             amount = 3),
                                             "OPLS-DA" = list(scaleC = data$info$dif$OPLSscaleC[[i]][j],
                                                              predI = 1,
                                                              orthoI = -1,
                                                              permI = ifelse(length(group) == 2,200,0),
                                                              log10L = data$info$dif$log10L[[i]][j],
                                                              amount = 3))
      
    }
  }
  
  if(all){
    group <- unique(data[["info"]][["sample"]][["class"]][[1]])
    compare$compare[["All"]] <- list(rawname = "All",
                                     group = unique(data[["info"]][["sample"]][["class"]][[1]]),
                                     "PCA" = list(scaleC = data$info$dif$PCAscaleC[[i]][j],
                                                  predI = -1,
                                                  orthoI = 0,
                                                  permI = 0,
                                                  log10L = data$info$dif$log10L[[i]][j],
                                                  amount = 3),
                                     "PLS-DA" = list(scaleC = data$info$dif$PLSscaleC[[i]][j],
                                                     predI = -1,
                                                     orthoI = 0,
                                                     permI = 0,
                                                     log10L = data$info$dif$log10L[[i]][j],
                                                     amount = 3),
                                     "OPLS-DA" = list(scaleC = data$info$dif$OPLSscaleC[[i]][j],
                                                      predI = 1,
                                                      orthoI = -1,
                                                      permI = 0,
                                                      log10L = data$info$dif$log10L[[i]][j],
                                                      amount = 3))
    group <- group[group!="QC"]
    compare$compare[["Allsample"]] <- list(rawname = "Allsample",
                                           group = group,
                                           "PCA" = list(scaleC = data$info$dif$PCAscaleC[[i]][j],
                                                        predI = -1,
                                                        orthoI = 0,
                                                        permI = 0,
                                                        log10L = data$info$dif$log10L[[i]][j],
                                                        amount = 3),
                                           "PLS-DA" = list(scaleC = data$info$dif$PLSscaleC[[i]][j],
                                                           predI = -1,
                                                           orthoI = 0,
                                                           permI = 200,
                                                           log10L = data$info$dif$log10L[[i]][j],
                                                           amount = 3),
                                           "OPLS-DA" = list(scaleC = data$info$dif$OPLSscaleC[[i]][j],
                                                            predI = 1,
                                                            orthoI = -1,
                                                            permI = 0,
                                                            log10L = data$info$dif$log10L[[i]][j],
                                                            amount = 3))
  }
  
  saveyaml(data = compare,filename = mulstatisticsyamlfile)
}