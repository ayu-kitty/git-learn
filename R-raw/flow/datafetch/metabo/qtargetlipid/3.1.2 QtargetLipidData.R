
#' @export
QTargetedLipidDataM <- function(data,
                                ...) {
  posdata <- QtargetLipidData(
    filename = data$info$datafrom$LC原始文件$posID,
    mode = "pos",
    grouping = data.frame(data$info$sample$LCMS$pos$name,
                          data$info$sample$samplename,
                          stringsAsFactors = F),
    weight = data$info$sample$weight)
  negdata <- QtargetLipidData(
    filename = data$info$datafrom$LC原始文件$negID,
    mode = "neg",
    grouping = data.frame(data$info$sample$LCMS$neg$name,
                          data$info$sample$samplename,
                          stringsAsFactors = F),
    weight = data$info$sample$weight)
  if (is.null(posdata) & is.null(negdata)) {
    stop("正负离子均为空，暂停分析")
  }
  
  mergedata <- rbind(posdata, negdata)
  mergedata <- mergedata[order(mergedata[, 3], -mergedata[, 5]), ]
  mergedata <- mergedata[!duplicated(mergedata$Lipid), ]
  data_final <- mutate(mergedata,
                       ID= 1:dim(mergedata)[1]) %>% 
    plyr::rename(c("Lipid"="Metabolites",
                   "Retention time"="Retention time (min)"))
  
  
  infodata <- getmetainfo(data = data_final[,"Metabolites",drop = F],idlist="Metabolites",
                          needlist = c("cid","CAS","HMDB","METLIN","Lipidmaps","KEGG","PubChem","ChEBI","smiles","InChIKey"))
  data_final <- merge(data_final,infodata,by = "Metabolites",sort = F)
  
  dataclass <- readdata(selectfile(path = "database/qualitative/QtargetLipid",file = "Lipid_match.xlsx"),sheet="QLipid")
  data_final <- merge(data_final,dataclass,by.x = "Metabolites",by.y = "Lipid")
  
  data_final_info <- data_final[,c("ID","m/z","Retention time (min)","Ion mode","Metabolites",
                                   "cid","CAS","HMDB","METLIN","Lipidmaps","KEGG","PubChem","ChEBI","smiles","InChIKey",
                                   "Class","Sub Class","Total Carbon","Total Unsaturation","Connection Method")]
  data_final_data <-  select(data_final,c("ID",data$info$sample$samplename))
  data_final <- merge(data_final_info,data_final_data,by = "ID")
  
  return(data_final)
}

#' @export
QtargetLipidData <- function(filename,
                             mode,
                             grouping,
                             weight) {
  if (any(is.na(weight))) {
    stop("请填写称重/体积列")
  } else if (any(weight == 0)) {
    stop("称重/体积列有0值")
  }
  
  if (is.na(filename)) {
    warning(paste0(mode, "模式下无定性文件位置信息"))
    return(NULL)
  } else if (!file.exists(filename)) {
    stop(paste0(mode, "模式下定性文件不存在"))
  }
  
  label <- readdata(filename = selectfile(path = "database/qualitative/QtargetLipid/",file = "拟靶向脂质数据库.xlsx"),
                    sheet = mode)
  
  label <- label[, c("Name", "TAG", "浓度", "Class")]
  data <- read.csv(file = filename,
                   sep = "\t", header = F,
                   row.names = 1, check.names = F,
                   stringsAsFactors = F, fileEncoding = "UTF-8", na.strings = "N/A")
  
  data <- data.frame(t(data), stringsAsFactors = F, check.names = F)
  data <- data[!is.na(data$`Identified metabolite`), ]
  data <- data[data$`Identified metabolite` != "Not annotated", ]
  data[is.na(data)] <- 0
  data$`Precursor m/z` <- as.numeric(data$`Precursor m/z`)
  data$`Product m/z` <- as.numeric(data$`Product m/z`)
  data$`Retention time` <- as.numeric(data$`Retention time`)
  data[, 5:dim(data)[2]] <- apply(data[, 5:dim(data)[2]], 2, as.numeric)
  data <- data[, -3]
  data <- merge(x = label, y = data, by.y = "Identified metabolite", by.x = "Name")
  
  
  # data <- data[,!grepl(pattern = "Blank",x = names(data))]
  data1 <- data[!is.na(data$TAG), ]
  data1 <- data1[!apply(X = (is.na(data1[, -1:-6])), MARGIN = 1, any), ]
  data1 <- data1[!apply(X = (data1[, -1:-6] == 0), MARGIN = 1, any), ]
  data2 <- data[is.na(data$TAG), ]
  
  if (any(grepl(pattern = "IS$", x = data2$Name))) {
    stop(paste("请确认", data2$Name[grepl(pattern = "IS$", x = data2$Name)],
               "是否在内标库中",
               collapse = ";"
    ))
  }
  
  for (i in 1:dim(data2)[1]) {
    class <- data2[i, "Class"]
    time <- data2[i, "Retention time"]
    label1 <- data1[data1$Class == class, ]
    
    if (dim(label1)[1] == 0) {
      label1 <- data1
      warning(paste0("在", mode, "中", data2[i, "Name"], "未找到内标"), immediate. = T)
    }
    
    if (dim(label1)[1] != 0) {
      if (any(label1$`Retention time` >= time)) {
        label1 <- label1[label1$`Retention time` >= time, ]
        label1 <- label1[which.min(label1$`Retention time`), ]
      } else {
        label1 <- label1[which.max(label1$`Retention time`), ]
      }
      
      # print(paste0(data2[i,"Name"],"使用",label1$Name,"作为内标"))
      
      data2[i, 7:dim(data2)[2]] <- data2[i, 7:dim(data2)[2]] / label1[1, 7:dim(label1)[2]] * label1$`浓度`
    } else {
      data2[i, "Precursor m/z"] <- NA
      warning(paste0("在", mode, "中", data2[i, "Name"], "未找到内标"), immediate. = T)
    }
  }
  
  data2[7:ncol(data2)] <- data2[7:ncol(data2)] * 0.2 / weight
  
  if (all(colnames(data2)[7:ncol(data2)] %in% grouping[, 1])) {
  } else {
    warning(paste0(mode, "模式下分析单多 或 定性少的样本："))
    print(colnames(data2)[7:ncol(data2)][colnames(data2)[7:ncol(data2)] %in% grouping[, 1]])
    stop(paste0(mode, "模式下样本信息不对应，暂停分析"))
  }
  
  if (all(grouping[, 1] %in% colnames(data2)[7:ncol(data2)])) {
  } else {
    warning(paste0(mode, "模式下分析单少 或 定性多的样本："))
    print(grouping[, 1][grouping[, 1] %in% colnames(data2)[7:ncol(data2)]])
    stop(paste0(mode, "模式下样本信息不对应，暂停分析"))
  }
  
  colnames(data2)[7:ncol(data2)] <- whto(a = grouping, b = colnames(data2)[7:ncol(data2)])
  
  data2 <- data2[, c(-2, -3, -4)]
  
  names(data2)[1:2] <- c("Lipid", "m/z")
  data2[, "Ion mode"] <- mode
  data2 <- data2[, c(2, dim(data2)[2], 1, 3, 4:(dim(data2)[2] - 1))]
  return(data2)
}

