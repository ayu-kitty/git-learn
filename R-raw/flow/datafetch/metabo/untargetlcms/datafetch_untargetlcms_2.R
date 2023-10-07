#!/opt/conda/bin/Rscript

#' LC非靶数据获取
#'
#' @param data obj
#' @param ... 见`datafetch_untargetlcms`
#'
#' @export
datafetch_untargetlcms_obj_2 <- function(data,
                                         ...){
  
  if(is.na(data$info$datafrom$LC原始文件$negM) & is.na(data$info$datafrom$LC原始文件$posM)){
    
    data$data$rawdata$data <- datafetch_untargetlcms_metapro(code.negID = data$info$datafrom$LC原始文件$negID,
                                                             code.posID = data$info$datafrom$LC原始文件$posID,
                                                             negname = data$info$sample$LCMS$neg$name,
                                                             posname = data$info$sample$LCMS$pos$name,
                                                             samplename =  data$info$sample$samplename,
                                                             ...)
    
    
  }else{
    data$data$rawdata$data <- datafetch_untargetlcms_2(code.negID = data$info$datafrom$LC原始文件$negID,
                                                       code.posID = data$info$datafrom$LC原始文件$posID,
                                                       code.negM = data$info$datafrom$LC原始文件$negM,
                                                       code.posM = data$info$datafrom$LC原始文件$posM,
                                                       negname = data$info$sample$LCMS$neg$name,
                                                       posname = data$info$sample$LCMS$pos$name,
                                                       samplename =  data$info$sample$samplename,
                                                       weight =  data$info$sample$weight,
                                                       ...)
  }
  
  return(data)
}

#' 获取QI定性定量文件
#' 
#' @param code.negID 负离子定性路径
#' @param code.posID 正离子定性路径 
#' @param code.negM 负离子M文件
#' @param code.posM 正离子M文件
#' @param negname 负离子M中样本名称 
#' @param posname 正离子M中样本名称 
#' @param samplename 样本分析名称
#' @param weight 样本称重
#' @param mode 处理模式，默认为"QI自带归一化模式"
#' @param deal mode为`raw`时,处理方式默认为"峰面积归一化"
#' @param neglab deal为`内标归一化`时,负离子内标名称
#' @param poslab deal为`内标归一化`时,正离子内标名称
#' @param data 定性数据
#' @param peptide 逻辑，是否保留肽段
#' @param meta 关注代谢物
#' @param score score筛选标准
#' @param fragscore 二级筛选标准
#' @param merge 逻辑值，是否仅保留未定性数据
#' @param ... 
#'
#' @export 
datafetch_untargetlcms_2 <- function(code.negID,
                                     code.posID,
                                     code.negM,
                                     code.posM,
                                     negname = NULL,
                                     posname = NULL,
                                     samplename = NULL,
                                     weight = NULL,
                                     mode = "normal",
                                     deal = "峰面积归一化",
                                     neglab = NULL,
                                     poslab = NULL,
                                     merge = F,
                                     ...) {
  
  
  negposID <- dataID_2(code.negID = code.negID,
                       code.posID = code.posID,
                       ...)
  
  negposM <- dataM(code.negM = code.negM,
                   code.posM = code.posM,
                   negname = negname,
                   posname = posname,
                   samplename = samplename,
                   weight = weight,
                   mode = mode, deal = deal, neglab = neglab, poslab = poslab)
  
  data <- MIDmerge(dataM = negposM,
                   dataID = negposID,
                   merge = merge)
  
  return(data)
}

#' 获取metapro定性定量文件
#' 
#' @param code.negID 负离子定性路径
#' @param code.posID 正离子定性路径 
#' @param negname 负离子M中样本名称 
#' @param posname 正离子M中样本名称 
#' @param samplename 样本分析名称
#' @param ... 
#'
#' @export 
datafetch_untargetlcms_metapro <- function(code.negID,
                                           code.posID,
                                           negname = NULL,
                                           posname = NULL,
                                           samplename = NULL,
                                           ...) {
  
  negdata <- getmetaprodata(filename = code.negID,
                            grouping = negname,
                            samplename= samplename,
                            ionmode = "neg")
  posdata <- getmetaprodata(filename = code.posID,
                            grouping = posname,
                            samplename= samplename,
                            ionmode = "pos")
  
  data <- rbind(negdata,posdata)
  data <- data[order(data[,dim(data)[2]],decreasing = T),]
  data <- data[!duplicated(data$cid),]
  data <- data[order(data$`Retention time (min)`),]
  data[,"ID"] <- 1:dim(data)[1]
  data <- data[,c(dim(data)[2],1:(dim(data)[2]-1))]
  
  return(data)
}

#' @export 
getmetaprodata <- function(filename,
                           grouping,
                           ionmode,
                           samplename){
  data2 <- readdata(filename = filename,sheet = "area")
  data2 <- data2[data2$`CheckStatus in set1` == "Success",]
  data2 <- data2[!duplicated(data2$name),]
  data2 <- data2[data2$`Mean RT in set1` > 0,]
  row.names(data2) <- data2$name
  data3 <- data2[,grepl(pattern = "\\[SAM\\]",x = colnames(data2)) | grepl(pattern = "\\[MIX\\]",x = colnames(data2))]
  colnames(data3) <- gsub(pattern = "\\[SAM\\]",replacement = "",x = colnames(data3))
  colnames(data3) <- gsub(pattern = "\\[MIX\\]",replacement = "",x = colnames(data3))
  data4 <- data3
  data4[is.na(data4)] <- 0
  
  if(!is.null(samplename)){
    if (all(grouping %in% colnames(data4))) {
      data4 <- data4[, grouping]
      colnames(data4) <- samplename
    } else {
      warning("样本量与样品信息不符合")
    } 
  }
  
  predata <- data2[,c("name","Mean RT in set1","ri","mainAdduct","mz")]
  colnames(predata) <- c("Compound ID","Retention time (min)","Average RI","Adducts","m/z")
  predata$Name <- gsub(pattern = "^[A-Za-z]*:",replacement = "",x = predata$`Compound ID`)
  predata <- getmetainfo(data = predata,idlist = "Name")
  row.names(predata) <- predata$`Compound ID`
  predata <- predata[!is.na(predata$cid),]
  predata[,"Ion mode"] <- ionmode
  
  predata <- predata[, c(
    "Retention time (min)","Average RI","Ion mode",
    "Name","cid", "HMDB","METLIN","Lipidmaps","KEGG","ChEBI","PubChem","CAS","smiles","InChIKey",
    "Super Class", "Class", "Sub Class",
    "Adducts", "Formula"
  )]
  colnames(predata)[colnames(predata) == "Name"] <- "Metabolites"
  predata <- merge(predata,data4,by = 0)[,-1]
}



