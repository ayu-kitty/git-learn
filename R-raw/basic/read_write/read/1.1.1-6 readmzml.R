#!/opt/conda/bin/Rscript

#' 读取mzml格式的文件 
#' 
#' @param filename 文件名 
#' @param ... 见`MSnbase::readSRMData`和`MSnbase::readMSData`
#' 
#' @export 
readmzml <- function(filename,
                     msLevel. = 1,
                     mode = "onDisk",
                     samplename = NULL,
                     ...){
  suppressMessages(library("MSnbase"))
  
  if(is.null(samplename)){
    samplename <- gsub(pattern = "\\.mzML$",replacement = "",x = basename(filename))
  }
  
  pdata <- AnnotatedDataFrame(data = data.frame(samplename = samplename,
                                                row.names = samplename,check.names = F))
  
  tryCatch({
    if(all(file.exists(filename))){
      
      data <-  MSnbase::readMSData(files = filename, 
                                   msLevel. = msLevel., 
                                   verbose = FALSE, 
                                   mode = mode,
                                   pdata = pdata,
                                   ...)
      
      if(dim(data@featureData@data)[1] == 0){
        data <-  MSnbase::readSRMData(files = filename,
                                      pdata = pdata)
      }
      
      return(data)
    }else{ 
      stop(paste0(filename,"文件不存在"))
    }
  })
}
