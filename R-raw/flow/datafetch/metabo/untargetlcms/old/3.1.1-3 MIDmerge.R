
#' @export
MIDmerge_old <- function(dataM, dataID, grouping, merge = T) {
  print("定量定量数据合并")
  negposM <- dataM
  negposID <- dataID
  negposMID <- merge(negposM, negposID,
    by.x = c("ID", "Ion mode"),
    by.y = c("ID", "Ion mode"),
    all.x = merge
  )
  negposMID2 <- subset(negposMID, select = c(ID, `m/z`, `Retention time (min)`, `Ion mode`, Metabolites:`Mass Error (ppm)`))
  negposMID3 <- negposMID[, names(negposMID) %in% grouping, drop = F]
  negposMID <- cbind(negposMID2, negposMID3)
  negposMID <- negposMID[order(-negposMID$Score), ]
  print("数据合并完毕")
  return(negposMID)
}
