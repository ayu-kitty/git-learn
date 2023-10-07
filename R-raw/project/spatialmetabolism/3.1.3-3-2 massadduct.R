#!/opt/conda/bin/Rscript

#' @export
adductMatchnew <- function(x,adductsdata,ppm=2, mDa=NULL) {
  # For each mass pair calculate the mass difference tolerance
  if (!is.null(mDa)) {
    Ad <- rep(mDa * 1e-3, times=length(x$A))
    Bd <- rep(mDa * 1e-3, times=length(x$B))
  } else {
    Ad <- x$A * ppm * 1e-6
    Bd <- x$B * ppm * 1e-6
  }
  # x$delta <- sqrt(Ad**2 + Bd**2) # Uncertainties add in quadrature
  x$delta <- apply(data.frame(Ad,Bd),1,min)
  
  indices <- vector()
  matches <- vector()
  
  adductsdata
  
  for ( i in 1:length(adductsdata$Name)) {
    mdata <- (x$A-adductsdata$Mass[i])*adductsdata$charge[i]/adductsdata$nM[i]
    
    for (j in 1:length(adductsdata$Name)) {
      if(i == j){
        next
      }
      mdata2 <- mdata*adductsdata$nM[j]/adductsdata$charge[j]+adductsdata$Mass[j]
      matchdiff <- abs(x$B - mdata2)
      idx <- which(matchdiff < x$delta)
      indices <- c(indices, idx)
      matches <- c(matches, rep(paste0(adductsdata$Name[i],";",adductsdata$Name[j]),length(idx)))
    }
  }
  
  output <- data.frame(A=x$A[indices],
                       B=x$B[indices],
                       diff=x$diff[indices],
                       delta=x$delta[indices],
                       matches=matches)
  
  class(output) <- c("massdiff","data.frame")
  row.names(output) <- NULL
  return(output)
}

#' @export
adductMatchnew2 <- function(x, add=mass2adduct::adducts, ppm=2, mDa=NULL) {
  # For each mass pair calculate the mass difference tolerance
  if (!is.null(mDa)) {
    Ad <- rep(mDa * 1e-3, times=length(x$A))
    Bd <- rep(mDa * 1e-3, times=length(x$B))
  } else {
    Ad <- x$A * ppm * 1e-6
    Bd <- x$B * ppm * 1e-6
  }
  # x$delta <- sqrt(Ad**2 + Bd**2) # Uncertainties add in quadrature
  x$delta <- apply(data.frame(Ad,Bd),1,min)
  
  indices <- vector()
  matches <- vector()
  mass <- vector()
  for (i in 1:length(add$mass)) {
    # Using a loop because number of adducts are few
    matchdiff <- abs(x$diff - add$mass[i])
    matchdiff2 <- x$diff - add$mass[i]
    idx <- which(matchdiff < x$delta)
    indices <- c(indices, idx)
    matches <- c(matches, rep(as.character(add$name[i]),length(idx)))
    mass <- c(mass,matchdiff2[idx]/x$delta[idx])
  }
  
  output <- data.frame(A=x$A[indices],
                       B=x$B[indices],
                       diff=x$diff[indices],
                       delta=x$delta[indices],
                       matches=matches,
                       mass=mass)
  # If original massdiff object contains correlation test results, include them too
  if (!is.null(x$Estimate)) {
    output$Estimate <- x$Estimate[indices]
  }
  if (!is.null(x$P.value)) {
    output$P.value <- x$P.value[indices]
  }
  if (!is.null(x$Significance)) {
    output$Significance <- x$Significance[indices]
  }
  
  class(output) <- c("massdiff","data.frame")
  row.names(output) <- NULL
  return(output)
}
