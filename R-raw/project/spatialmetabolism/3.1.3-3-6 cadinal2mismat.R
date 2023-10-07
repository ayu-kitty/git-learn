
#' @export
cardinal2msimat2 <- function(d) UseMethod("cardinal2msimat2")

#' @export
cardinal2msimat2.MSContinuousImagingExperiment <- function(d) {
  # Check if peaks have been binned, otherwise the keys will not correspond to
  # peaks
  # if ("peakBin" %in% names(d@processing)) {
  if (TRUE) {
    # Check that data are indeed in sparse matrix format
    if ("matrix" %in% class(spectra(d))) {
      # Initialize vectors for containing data
      peaks <- mz(d) # Binned mass values
      spots <- 1:dim(d)['Pixels'] # Pixel name values
      data <- t(spectra(d))
      peakintensities <- colSums(data)
      # Store as list and declare class as msimat
      out <- list(mat=as.matrix(data),
                  peaks=as.numeric(peaks),
                  spots=spots,
                  peakintensities=peakintensities)
      class(out) <- "msimat"
      # Return msimat object
      return(out)
    } else {
      stop("Matrix not stored internally as sparse_matc, which should be the case for MSProcessedImagingExperiment")
    }
  } else {
    stop("Only applicable to binned data. Please apply peakBin() processing first.")
  }
}

#' @export
cardinal2msimat2.MSProcessedImagingExperiment <- function(d) {
  
  # Check if peaks have been binned, otherwise the keys will not correspond to
  # peaks
  # if ("peakBin" %in% names(d@processing)) {
  if (TRUE) {
    # Check that data are indeed in sparse matrix format
    if ("sparse_matc" %in% class(spectra(d))) {
      # Initialize vectors for containing data
      peaks <- spectra(d)@keys # Binned mass values
      spots <- 1:dim(d)['Pixels'] # Pixel name values
      cols <- vector()
      rows <- vector()
      vals <- vector()
      # For each spot, parse the relevant vector
      for (i in 1:length(spots)) {
        vals <- c(vals, spectra(d)@data$values[[i]])
        rows <- c(rows, rep(i-1, length(spectra(d)@data$values[[i]])))
        peaksvec <- spectra(d)@data$keys[[i]]
        peaksidx <- sapply(peaksvec, function(p) which (peaks==p) - 1)
        cols <- c(cols, peaksidx) # -1 to convert to 0-based indexing
      }
    } else {
      stop("Matrix not stored internally as sparse_matc, which should be the case for MSProcessedImagingExperiment")
    }
  }
  else {
    stop("Only applicable to binned data. Please apply peakBin() processing first.")
  }
  
  # Convert to sparseMatrix object
  tsm <- Matrix::sparseMatrix(i=rows,
                              j=cols,
                              x=vals,
                              index1=FALSE,
                              giveCsparse=FALSE,
                              check=TRUE)
  peakintensities <- Matrix::colSums(tsm)
  # Store as list and declare class as msimat
  data <- list(mat=tsm,
               peaks=peaks,
               spots=spots,
               peakintensities=peakintensities)
  class(data) <- "msimat"
  # Return msimat object
  return(data)
}