
rbind_lm <- function(...){
  data <- list(...)
  allcolname <- NULL
  for ( i in 1:length(data)) {
    allcolname <- c(allcolname,colnames(data[[i]]))
  }
  allcolname <- unique(allcolname)
  
  for ( i in 1:length(data)) {
    for ( j in 1:length(allcolname)) {
      if (!(allcolname[j] %in% colnames(data[[i]]))) {
        newdata <- data[[i]]
        newdata[,allcolname[j]] <- NA
        data[[i]] <- newdata
      }
    }
  }
  
  data <- do.call(rbind,data)
}

QtargetLipid <- data.frame(type = "拟靶向脂质",
                           S3 = "QtargetLipid",
                           row.names = "拟靶向脂质")

QtargetMic <- data.frame(type = "拟靶向肠道菌群",
                         S3 = "QtargetMic",
                         row.names = "拟靶向肠道菌群")

SpaceM <- data.frame(type = "空间代谢组",
                     S3 = "SpaceM",
                     row.names = "空间代谢组")

#' @export
projectinfo <- rbind_lm(QtargetLipid,
                        QtargetMic,
                        SpaceM)
