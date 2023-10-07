#!/opt/conda/bin/Rscript

#' M和ID数据合并
#' 
#' @param dataM M数据
#' @param dataID ID数据
#' @param samplename 样本分析名 
#' @param merge 逻辑值，是否仅保留未定性数据
#'
#' @export
MIDmerge <- function(dataM, 
                     dataID, 
                     merge = T) {
  print("定量定量数据合并")
  negposM <- dataM
  negposID <- dataID
  negposMID <- merge(negposM, negposID,
                     by.x = c("ID", "Ion mode"),
                     by.y = c("ID", "Ion mode"),
                     all.x = merge)
  negposMID2 <- subset(negposMID, select = c(ID, `m/z`, `Retention time (min)`,
                                             `Ion mode`, Metabolites:`Mass Error (ppm)`))
  negposMID3 <- negposMID[, !(colnames(negposMID) %in% colnames(negposMID2)), drop = F]
  negposMID <- cbind(negposMID2, negposMID3)
  negposMID <- negposMID[order(-negposMID$Score), ]
  negposMID2 <- negposMID
  negposMID$Metabolites <- stringr::str_to_title(negposMID$Metabolites)
  negposMID$Metabolites[which((grepl("^[A-Za-z]+.*\\(.*:.*\\)$", negposMID$Metabolites)))] <- negposMID2$Metabolites[which((grepl("^[A-Za-z]+.*\\(.*:.*\\)$", negposMID$Metabolites)))]
  print("数据合并完毕")
  return(negposMID)
}
