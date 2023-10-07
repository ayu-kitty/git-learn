
#' @export
smooth.image.gaussian <- function(x, window=3, ...) {
  if ( all(is.na(x)) ) return(x)
  r <- floor(window / 2)
  sd <- window / 4
  x.new <- .Call("C_gaussianFilter", x, r, sd, PACKAGE="Cardinal")
  x.new <- max(x, na.rm=TRUE) * x.new / max(x.new, na.rm=TRUE)
  x.new
}

#' @export
smooth.image.adaptive <- function(x, window=3, ...) {
  if ( all(is.na(x)) ) return(x)
  r <- floor(window / 2)
  sd <- window / 4
  x.new <- .Call("C_bilateralFilter", x, r, sd, PACKAGE="Cardinal")
  x.new <- max(x, na.rm=TRUE) * x.new / max(x.new, na.rm=TRUE)
  x.new
}
