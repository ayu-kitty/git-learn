
#' 空代数据条纹处理
#'
#' @param mse 
#'
#' @export
imzmlstripdeal <- function(mse,saverawdata = F){
  loadpython()
  stripdeal <- import("lmbio.spatialmetabolism.fixstrip")
  
  spectradata <- as.matrix(iData(mse, "intensity"))
  spectradata <- as.data.frame(spectradata)
  colnames(spectradata) <- paste(coord(mse)$x, coord(mse)$y,
                                 sep = "-")
  row.names(spectradata) <- format(mz(mse), nsmall = 5, trim = T)
  
  if(saverawdata){
    savetxt(data = spectradata,filename = "strip.txt",row.names = TRUE)
  }
  
  spectradata2 <- stripdeal$fixstrip_data(data = spectradata)
  spectradata2[spectradata == 0] <- 0
  spectradata2 <- as.matrix(spectradata2)
  spectra(mse)[,] <- spectradata2
  return(mse)
}
