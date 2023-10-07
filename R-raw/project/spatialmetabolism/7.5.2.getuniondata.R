#' getuniondata
#'
#' 将空代和空转数据进行整理
#'
#' @param samplename 样本名称
#' @param mode 模式
#' @param metadatapath 空代数据路径
#' @param metaqualitative 空代定性数据路径
#' @param transdatapath 空转数据路径
#' @param savepath 结果保存路径
#'
#' @export
getuniondata <- function(samplename,
                         mode = "neg",
                         metadatapath = "./sample/union/data/",
                         metaqualitative = "./sample/qualitative/Qualitative.xlsx",
                         transdatapath = "./sample/union/transdata/",
                         savepath = "./sample/union/uniondata/"){
  
  transdata <- read.csv(file = paste0(transdatapath,samplename,"/","counts.csv"), header = T, check.names = F)
  tissuedata <- read.csv(file = paste0(transdatapath,samplename,"/","tissue_positions_list.csv"),header = F)[,1:4]
  colnames(tissuedata) <- c("transname","sample","y","x")
  tissuedata[,"metaname"] <- paste(samplename,mode,tissuedata[,"x"],tissuedata[,"y"],sep = "-")
  tissuedata <- tissuedata[tissuedata[,"sample"] == 1,]

  metadata <- read.table(file = paste0(metadatapath,samplename,"-",mode,".txt"),
                         header = T,check.names = F,stringsAsFactors = F)
  
  # 统一barcode标识，LN3.AAACAAGTATCTCCCA
  tissuedata[,"transname"] <- paste0(samplename, ".", gsub(pattern = "-[0-9]", replacement = "", tissuedata[,"transname"]))
  tissuedata <- tissuedata[tissuedata[,"metaname"] %in% colnames(metadata)[-1],]
  tissuedata <- tissuedata[tissuedata[,"transname"] %in% colnames(transdata)[-1],]

  # 提取转录数据，并整理
  transdata2 <- transdata[,c("Gene",tissuedata[,"transname"])]
  transdata2 <- data.frame("Accession" = transdata[,c("Gene")], transdata2, check.names = F)
  colnames(transdata2)[2] <- "Gene Name"

  # 提取代谢数据，并整理
  metadata2 <- metadata[,c("mz",tissuedata[,"metaname"])]
  qualitativedata <- lmbio::readdata(filename = metaqualitative,sheet = mode)
  metadata2 <- cbind(qualitativedata,metadata2[,-1])

  if (!dir.exists(paste(savepath,samplename,mode,sep = "/"))) {
    dir.create(path = paste(savepath,samplename,mode,sep = "/"), recursive = T)
  }

  write.table(x = tissuedata,
              file = paste(savepath,samplename,mode,"samplename.txt",sep = "/"),
              sep = "\t", row.names = F)
  write.table(x = transdata2,
              file = paste(savepath,samplename,mode,"transdata.txt",sep = "/"),
              sep = "\t", row.names = F)
  write.table(x = metadata2,
              file = paste(savepath,samplename,mode,"metadata.txt",sep = "/"),
              sep = "\t", row.names = F)
}

