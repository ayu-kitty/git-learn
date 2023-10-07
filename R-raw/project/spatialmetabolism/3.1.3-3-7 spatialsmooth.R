projectToRaster2 <- function(x, y, values, dim = NULL, res = NULL) {
  xrange <- range(x, na.rm=TRUE)
  yrange <- range(y, na.rm=TRUE)
  if ( is.null(dim) ) {
    asp <- diff(yrange) / diff(xrange)
    nrow <- as.integer(sqrt(length(values) / asp))
    ncol <- as.integer(asp * nrow)
  } else {
    nrow <- dim[1]
    ncol <- dim[2]
  }
  if ( is.null(res) || is.na(res[1]) ) {
    rx <- (xrange[2] - xrange[1]) / (nrow - 1)
  } else {
    rx <- res[1]
  }
  if ( is.null(res) || is.na(res[2]) ) {
    ry <- (yrange[2] - yrange[1]) / (ncol - 1)
  } else {
    ry <- res[2]
  }
  rx <- ifelse(is.finite(rx), rx, 1)
  ry <- ifelse(is.finite(ry), ry, 1)
  rows <- as.integer(round((x - xrange[1]) / rx))
  cols <- as.integer(round((y - yrange[1]) / ry))
  init <- as.vector(NA, mode=typeof(values))
  rs <- matrix(init, nrow=nrow, ncol=ncol)
  idx <- rows + cols * nrow + 1L
  valid <- is.finite(idx)
  values <- values[valid]
  idx <- idx[valid]
  idx[idx > length(rs)] <- length(rs)
  rs[idx] <- values
  rs
}

spatialsmoothdata <- function(values,x,y,window = 2,smooth.image = "gaussian"){
  x <- x-min(x)+1
  y <- y-min(y)+1
  tproj <- projectToRaster2(x = x, y = y, values = values, dim = c(max(x),max(y)), res = NULL)
  tproj <- structure(tproj, range=c(NA,NA), resolution = 1)
  if(smooth.image == "gaussian"){
    tproj <- smooth.image.gaussian(tproj,window = window)
  }else if(smooth.image == "adaptive"){
    tproj <- smooth.image.adaptive(tproj,window = window)
  }else{
    stop(paste0("smooth.image不支持",smooth.image,"方法"))
  }
  for ( i in 1:length(values)) {
    if(!is.na(tproj[x[i],y[i]])){
      values[i] <- tproj[x[i],y[i]]
    }
  }
  return(values)
}

#' @export
spatialsmooth <- function(.object,window = 2,smooth.image = "gaussian"){
  
  x <- coord(.object)$x
  y <- coord(.object)$y
  
  if(smooth.image == "none"){
    spectradata2 <- spectra(.object)
  }else{
    spectradata2 <- featureApply(.object = .object,.fun = spatialsmoothdata,x = x,y = y,
                                 window = window,smooth.image = smooth.image)
    spectradata2 <- t(spectradata2)
    spectradata2 <- as.matrix(spectradata2)
  }
  
  return(spectradata2)
}
