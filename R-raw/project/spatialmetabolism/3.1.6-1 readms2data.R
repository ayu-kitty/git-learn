#!/opt/conda/bin/Rscript

#提取数据
readms2data <- function(path,precursorMZ){
  suppressMessages(library("readMzXmlData"))
  ms2data <- dir(path,full.names = TRUE,recursive = FALSE)
  spec <- readMzXmlDir(file.path(path))
  print("读取全部RAW文件完成...")
  spec_ion <- list()
  #从所有spectrum中提取precursorMz在范围内的spectrum
  for(i in c(1:length(spec))){
    if(abs(spec[[i]][["metaData"]][["precursorMz"]]-precursorMZ)/precursorMZ*10^6<20){
      spec_ion[i] <-spec[i]}
  }
  spec_ion <- spec_ion[!sapply(spec_ion,is.null)]
  #将第112次扫描去掉
  if(length(spec_ion)==5600){
    for(i in c(1:50)){
      spec_ion[[i*112]] <- NA
    }
  }
  spec_ion <- spec_ion[!is.na(spec_ion)]
  
  mass <- numeric()
  intensity <- numeric()
  samples <- numeric()
  #提取该离子所有scan的的mz & intensity
  for(i in 1:length(spec_ion)){
    sample1 <- rep(i,length(spec_ion[[i]][["spectrum"]][["mass"]]))
    samples <- c(samples,sample1)
    mass <- c(mass,spec_ion[[i]][["spectrum"]][["mass"]])
    intensity <- c(intensity,spec_ion[[i]][["spectrum"]][["intensity"]])
  }
  
  spec_1 <- cbind(samples,mass,intensity)
  spec_1 <- data.frame(spec_1)
  return(spec_1)
}

#峰校正
binpeak <- function(data,tolerance=10e-6){
  samples <- data[ ,1]
  mass <- data[ ,2]
  intensities <-data[ ,3]
  ## sort values by mass
  s <- sort.int(mass, index.return=TRUE)
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  
  n <- length(mass)
  
  d <- diff(mass)
  nBoundaries <- max(100L, floor(3L * log(n)))
  boundary <- list(left=double(nBoundaries), right=double(nBoundaries))
  
  currentBoundary <- 1L
  boundary$left[currentBoundary] <- 1L
  boundary$right[currentBoundary] <- n
  
  while (currentBoundary > 0L) {
    ## find largest gap
    left <- boundary$left[currentBoundary]
    right <- boundary$right[currentBoundary]
    currentBoundary <- currentBoundary - 1L
    gaps <- d[left:(right-1L)]
    
    gapIdx <- which.max(gaps) + left - 1L
    
    ## left side
    
    l_mass=mass[left:gapIdx]
    l_intensities=intensities[left:gapIdx]
    l_samples=samples[left:gapIdx]
    l_meanMass <- mean(l_mass)
    ## further splitting needed?
    if (any(abs(l_mass - l_meanMass) / l_meanMass > tolerance)) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- left
      boundary$right[currentBoundary] <- gapIdx
    } else {
      mass[left:gapIdx] <-  l_meanMass
      #intensities[left:gapIdx] <-max(intensities[left:gapIdx])
    }
    
    ## right side
    r_mass=mass[(gapIdx + 1L):right]
    r_intensities=intensities[(gapIdx + 1L):right]
    r_samples=samples[(gapIdx + 1L):right]
    r_meanMass <- mean(r_mass)
    ## further splitting needed?
    if (any(abs(r_mass - r_meanMass) / r_meanMass > tolerance)) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- gapIdx + 1L
      boundary$right[currentBoundary] <- right
    } else {
      mass[(gapIdx + 1L):right] <- r_meanMass
      #intensities[(gapIdx + 1L):right] <- max(intensities[(gapIdx + 1L):right])
    }
  }
  ## stack size have to be increased?
  ## (should rarely happen because recursion deep is mostly < 20)
  #if (currentBoundary == nBoundaries) {
  #nBoundaries <- floor(nBoundaries * 1.5)
  #boundary$left <- c(boundary$left,
  #double(nBoundaries - currentBoundary))
  #boundary$right <- c(boundary$right,
  # double(nBoundaries - currentBoundary))
  
  data_bin <- data.frame(samples,mass,intensities)
  data_bin <- data_bin[order(data_bin$samples),]
  return(data_bin)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  filepath = "F:/同位素处理/code2/二级_分步/PRM-mzxml/"
  precursorMZ = 146.1649
  # precursorMZ = 203.2227
  # precursorMZ = 307.0433
  # precursorMZ = 798.5390
  # precursorMZ = 826.5705
  outfile = paste0(precursorMZ,"_bin.csv")
  
  spec_1 <- readms2data(path=filepath,precursorMZ = precursorMZ)
  # write.csv(paste0(precursorMZ,".csv"),file=outfile,row.names = FALSE)
  data_bin <- binpeak(spec_1)
  print("峰校正完成...")
  write.csv(data_bin,file=outfile,row.names = FALSE)
}
