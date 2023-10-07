#' 补全空代数据
#'
#' @export
imzmlfill <- function(mse){
  
  mse <- Cardinal::pull(mse, as.matrix = TRUE)
  coorddata <- as.data.frame(coord(mse))
  mse <- mse[,!duplicated(coorddata)]
  coorddata <- as.data.frame(coord(mse))
  coorddata2 <- expand.grid(x=unique(coorddata$x), y=unique(coorddata$y))
  coorddata2 <- coorddata2[!(paste(coorddata2$x,coorddata2$y) %in% paste(coorddata$x,coorddata$y)),]
  
  if(dim(coorddata2)[1] > 0){
    pixelData <- as.data.frame(pData(mse))
    if("X3DPositionX" %in% colnames(pixelData)){
      pixelData2 <- PositionDataFrame(coord = coorddata2,run = pixelData[1,"run"],
                                      "3DPositionX" = coorddata2$x,
                                      "3DPositionY" = coorddata2$y,
                                      pixelData[1,!(colnames(pixelData) %in% c("run","x","y","X3DPositionX","X3DPositionY")),drop = F],
                                      check.names = F)
    }else{
      pixelData2 <- PositionDataFrame(coord = coorddata2,
                                      run = pixelData[1,"run"],
                                      pixelData[1,!(colnames(pixelData) %in% c("run","x","y"))],
                                      check.names = F)
    }
    
    bindmse <- MSImagingExperiment(imageData = matrix(rep(0,dim(coorddata2)[1]*length(mz(mse))),
                                                      nrow = length(mz(mse)),ncol = dim(coorddata2)[1]),
                                   featureData = MassDataFrame(mz = mz(mse)),
                                   pixelData = pixelData2,
                                   centroided = T)
    
    mse <- BiocGenerics::cbind(mse, bindmse)
  }
  
  mse <- mse[, order(coord(mse)$y, coord(mse)$x)]
  
  return(mse)
}
