#!/opt/conda/bin/Rscript

# getmatchformula
#' @export
getmatchformula <- function(mzdata,
                            formuladata,
                            mode,
                            ppm = 10) {
  
  mz <- mzdata
  acumeta <- formuladata
  acumeta <- acumeta[acumeta$mode == mode,]
  
  matchmz2 <- data.frame()
  # 匹配
  for (i in 1:length(mz)) {
    minmz <- mz[i] - ppm / 1000000 * mz[i]
    maxmz <- mz[i] + ppm / 1000000 * mz[i]
    matchmz <- acumeta[((acumeta$mz <= maxmz) & (acumeta$mz >= minmz)), ]
    matchmz[,"ppm"] <- abs((matchmz$mz-mz[i])/matchmz$mz)*1000000
   
    if (dim(matchmz)[1] > 0) {
      matchmz <- matchmz[order(matchmz$mz),]
      matchmz[,"mz"] <- mz[i]
      matchmz2 <- rbind(matchmz2,matchmz)
    }
  }
  
  return(matchmz2)
}