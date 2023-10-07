#!/opt/conda/bin/Rscript

#' @export
getareadatainreport <- function(qualitativefile = "../sample/qualitative/Qualitative.xlsx",
                                areadatapath = "../sample/area/data/",
                                savepath = "./3.选区数据/数据矩阵",
                                moderange = c("neg", "pos")){
  
  if(!dir.exists(areadatapath)){
    warning(paste0(areadatapath,"目录为空，跳过运行"),immediate. = T)
    return()
  }
  
  for (mode in moderange) {
    filesheet <- getsheetname(qualitativefile)
    if(!(mode %in% filesheet)){
      next
    }
    
    qdata <- readdata(qualitativefile,
                      sheet = mode)
    
    filename <- list.files(path = areadatapath, 
                           pattern = paste0("-", mode, "-pixel_level.txt"),
                           full.names = T)
    
    for (i in seq_len(length(filename))) {
      data <- read.table(file = filename[i],sep = "\t",header = T,check.names = F)
      mzdata <- cbind(qdata, data[, -1, drop = F])
      
      savexls(data = mzdata,
              filename = paste0(savepath,"/",gsub(pattern = ".txt",replacement = paste0("-",mode,".xls"),x = basename(filename[i]))))
    }
    
    
    filename <- list.files(path = areadatapath,
                           pattern = paste0("-", mode, "-sample_level.txt"),
                           full.names = T)
    
    data <- read.table(file = filename[1],sep = "\t",header = T,check.names = F)
    
    if (length(filename) > 1) {
      for (i in 2:length(filename)) {
        data1 <- read.table(file = filename[i],sep = "\t",header = T,check.names = F)
        data <- cbind(data, data1[, -1, drop = F])
      }
    }
    
    mzdata <- cbind(qdata, data[, -1, drop = F])
    
    savexls(data = mzdata,
            filename = paste0(savepath,"/Allsample-", mode, "-sample_level.xls"))
  }
}

