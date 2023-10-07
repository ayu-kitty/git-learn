# Intestine Pseudo Targeting

#' @export
micfetch <- function(data) {
  data1 <- data
  GCdata <- IntestineGCfileDeal(
    GCFile = data1$info$datafrom$GC原始文件,
    name = data.frame(data1$info$sample$GCMS$name, data1$info$sample$samplename)
  )
  LCdata_neg <- IntestineLCfileDeal(
    LCFile = data1$info$datafrom$LC原始文件$negID,
    mode = "neg",
    name = data.frame(data1$info$sample$LCMS$neg$name, data1$info$sample$samplename)
  )
  LCdata_pos <- IntestineLCfileDeal(
    LCFile = data1$info$datafrom$LC原始文件$posID,
    mode = "pos",
    name = data.frame(data1$info$sample$LCMS$pos$name, data1$info$sample$samplename)
  )
  LCdata <- rbind(LCdata_neg, LCdata_pos)
  rawdata <- rbind(LCdata, GCdata)
  rawdata[, "name"] <- tolower(rawdata$`Component Name`)

  METAdata1 <- getdatabase()
  METAdata2 <- getdatabase(database = "info.db")[, c("id", "super_class", "class", "sub_class", "kegg")]
  METAdata <- merge(METAdata1, METAdata2)
  names(METAdata) <- c("Compound ID", "name", "Super Class", "Class", "Sub Class", "kegg")
  rawdata <- merge(METAdata, rawdata, by.x = "name", by.y = "name", all.y = T)

  rawdata$name <- rawdata$`Component Name`
  rawdata <- rawdata[, names(rawdata) != "Component Name"]
  names(rawdata)[1] <- "Metabolites"

  rawdata <- rawdata[order(apply(rawdata[, -1:-8], 1, mean), decreasing = T), ]
  rawdata <- rawdata[!duplicated(rawdata$Metabolites), ]

  data1$data$rawdata$data <- rawdata

  return(data1)
}


#' @export
IntestineGCfileDeal <- function(GCFile,
                                name) {

  # GClabel <- RMETA2::readxlsx(name = "数据库.xlsx",sheet = "GC-内标")
  # GCdatabase <- RMETA2::readxlsx(name = "数据库.xlsx",sheet = "GC-数据库")

  data <- readdata(GCFile)
  data <- reshape2::dcast(data = data, Compound + RT ~ Filename, value.var = "Area")
  data <- tidyr::separate(data = data, col = "RT", into = "RT", sep = ",")
  data <- data[data$RT != 0, ]
  data <- data[, !grepl(pattern = "[Bb][Ll][Aa][Nn][Kk]", x = names(data))]
  # data <- data[,!grepl(pattern = "0603blk",x = names(data))]
  # data <- data[,!grepl(pattern = "0603mix",x = names(data))]
  data[data == "N/F"] <- 0
  data[data == "N/A"] <- 0
  data[is.na(data)] <- 0
  data[, -1] <- apply(data[, -1], 2, as.numeric)

  # datalabel <- data[data$Compound %in% GClabel$Name,]
  # datalabel <- datalabel[order(RMETA2::whto(a = GClabel[,c("Name","Level")],b = datalabel$Compound),decreasing = T),]
  # datalabel <- tidyr::separate(data = datalabel,col = "Compound",into = "Compound",sep = " _1")
  # datalabel <- tidyr::separate(data = datalabel,col = "Compound",into = "Compound",sep = " _2")
  # datalabel <- tidyr::separate(data = datalabel,col = "Compound",into = "Compound",sep = " _3")
  # datalabel <- tidyr::separate(data = datalabel,col = "Compound",into = "Compound",sep = " _4")
  # datalabel <- tidyr::separate(data = datalabel,col = "Compound",into = "Compound",sep = "_1")
  # datalabel <- tidyr::separate(data = datalabel,col = "Compound",into = "Compound",sep = "_2")
  # datalabel <- tidyr::separate(data = datalabel,col = "Compound",into = "Compound",sep = "_3")
  # datalabel <- tidyr::separate(data = datalabel,col = "Compound",into = "Compound",sep = "_4")
  # datalabel <- datalabel[!duplicated(datalabel$Compound),]
  # datalabel <- datalabel[order(datalabel$RT),]
  #
  # if(any(datalabel[,3:dim(datalabel)[2]] == 0)){
  #
  #   datalabel <- datalabel[apply(X = datalabel[,3:dim(datalabel)[2]],MARGIN = 1,FUN = function(x){all(x!=0)}),]
  #   warning("内标中有0值")
  #
  # }
  #
  # data <- data[data$Compound %in% GCdatabase$Name,]
  # data <- data[order(RMETA2::whto(a = GCdatabase[,c("Name","Level")],b = data$Compound),decreasing = T),]

  data[, "Compound"] <- gsub(pattern = "-1$", replacement = "", x = data[, "Compound"])
  data[, "Compound"] <- gsub(pattern = "-2$", replacement = "", x = data[, "Compound"])
  data <- tidyr::separate(data = data, col = "Compound", into = "Compound", sep = " _1")
  data <- tidyr::separate(data = data, col = "Compound", into = "Compound", sep = " _2")
  data <- tidyr::separate(data = data, col = "Compound", into = "Compound", sep = " _3")
  data <- tidyr::separate(data = data, col = "Compound", into = "Compound", sep = " _4")
  data <- tidyr::separate(data = data, col = "Compound", into = "Compound", sep = "_1")
  data <- tidyr::separate(data = data, col = "Compound", into = "Compound", sep = "_2")
  data <- tidyr::separate(data = data, col = "Compound", into = "Compound", sep = "_3")
  data <- tidyr::separate(data = data, col = "Compound", into = "Compound", sep = "_4")
  data <- data[order(apply(data[, -1:-2], 1, mean), decreasing = T), ]
  data <- data[apply(data[, -1:-2], 1, mean) != 0, ]
  data <- data[!duplicated(data$Compound), ]
  data <- data[order(data$RT), ]


  # for ( i in 1:dim(data)[1]) {
  #
  #   label2 <- datalabel[datalabel$RT > data$RT[i],]
  #
  #   if(dim(label2)[1] == 0){
  #     label2 <- datalabel[which.max(datalabel$RT),]
  #   }else{
  #     label2 <- label2[which.min(label2$RT),]
  #   }
  #
  #   data[i,3:dim(data)[2]] <- data[i,3:dim(data)[2]] / label2[,3:dim(data)[2]]
  #
  # }

  i <- 3
  while (i <= dim(data)[2]) {
    data[, i] <- data[, i] / sum(data[, i]) * 10000
    i <- i + 1
  }

  colnames(data)[1:2] <- c("Component Name", "Retention time")
  data[, "mode"] <- "GC"
  data <- data[, c(1, 2, dim(data)[2], 4:dim(data)[2] - 1)]
  data <- data[, c(T, T, T, (colnames(data)[-1:-3] %in% name[, 1]))]
  colnames(data)[-1:-3] <- whto(name, colnames(data)[-1:-3])

  return(data)
}

#' @export
IntestineLCfileDeal <- function(LCFile,
                                mode,
                                name) {
  LClabel <- readdata(name = "数据库.xlsx", sheet = paste0("LC-", mode, "-内标"))
  LCdatabase <- readdata(name = "数据库.xlsx", sheet = paste0("LC-", mode, "-数据库"))

  data <- read.csv(file = LCFile, sep = "\t", header = T, check.names = F, stringsAsFactors = F, fileEncoding = "UTF-8")

  if (names(data)[1] == "Index") {
    data <- data[, c("Sample Name", "Component Name", "Area", "Retention Time")]
    data <- data[data[, 3] != "N/A", ]
    data[is.na(data)] <- 0
    data[, 3] <- as.numeric(data[, 3])
    data[, 4] <- as.numeric(data[, 4])
    data1 <- reshape2::dcast(data = data, formula = `Component Name` ~ `Sample Name`, value.var = "Area", fill = 0)
    data2 <- reshape2::dcast(data = data, formula = `Component Name` ~ "Retention Time", fun.aggregate = mean, value.var = "Retention Time")
    data <- merge(data2, data1)
    names(data)[2] <- "Retention time"
  } else {
    data <- read.csv(file = LCFile, sep = "\t", header = F, row.names = 1, check.names = F, stringsAsFactors = F, fileEncoding = "UTF-8", na.strings = "N/A")
    data <- data.frame(t(data), stringsAsFactors = F, check.names = F)
    data <- data[!is.na(data$`Identified metabolite`), ]
    data <- data[data$`Identified metabolite` != "Not annotated", ]
    data[is.na(data)] <- 0
    data$`Precursor m/z` <- as.numeric(data$`Precursor m/z`)
    data$`Product m/z` <- as.numeric(data$`Product m/z`)
    data$`Retention time` <- as.numeric(data$`Retention time`)
    data[, 5:dim(data)[2]] <- apply(data[, 5:dim(data)[2]], 2, as.numeric)
    data <- data[, -2]
    data <- data[, -3]
    names(data)[1:2] <- c("Component Name","Retention time")
    data <- data[, c(1, 2, 3:dim(data)[2])]
  }

  data <- data[, !grepl(pattern = "[Bb][Ll][Aa][Nn][Kk]", x = names(data))]
  # data <- data[,!grepl(pattern = "blk-inject",x = names(data))]
  # data <- data[,!grepl(pattern = "blk-pre",x = names(data))]
  # data <- data[,!grepl(pattern = "STD148mix",x = names(data))]

  datalabel <- data[data$`Component Name` %in% LClabel$ID, ]
  # datalabel <- datalabel[order(RMETA2::whto(a = LClabel[,c("ID","Level")],b = datalabel$`Component Name`),decreasing = T),]
  datalabel <- tidyr::separate(data = datalabel, col = "Component Name", into = "Component Name", sep = "-1")
  datalabel <- tidyr::separate(data = datalabel, col = "Component Name", into = "Component Name", sep = "-2")
  datalabel <- tidyr::separate(data = datalabel, col = "Component Name", into = "Component Name", sep = "-3")
  datalabel <- tidyr::separate(data = datalabel, col = "Component Name", into = "Component Name", sep = "-4")
  datalabel <- datalabel[order(apply(datalabel[, -1:-2], 1, mean), decreasing = T), ]
  datalabel <- datalabel[apply(datalabel[, -1:-2], 1, mean) != 0, ]
  datalabel <- datalabel[!duplicated(datalabel$`Component Name`), ]
  datalabel <- datalabel[order(datalabel$`Retention time`), ]

  if (any(datalabel[, 3:dim(datalabel)[2]] == 0)) {

    datalabel <- datalabel[apply(X = datalabel[, 3:dim(datalabel)[2]], MARGIN = 1, FUN = function(x) {
      all(x != 0)
    }), ]
    warning("内标中有0值")
  }


  data <- data[data$`Component Name` %in% LCdatabase$ID, ]
  # data <- data[order(RMETA2::whto(a = LCdatabase[,c("ID","Level")],b = data$`Component Name`),decreasing = T),]
  data <- tidyr::separate(data = data, col = "Component Name", into = "Component Name", sep = "-1")
  data <- tidyr::separate(data = data, col = "Component Name", into = "Component Name", sep = "-2")
  data <- tidyr::separate(data = data, col = "Component Name", into = "Component Name", sep = "-3")
  data <- tidyr::separate(data = data, col = "Component Name", into = "Component Name", sep = "-4")
  data <- data[order(apply(data[, -1:-2], 1, mean), decreasing = T), ]
  data <- data[apply(data[, -1:-2], 1, mean) != 0, ]
  data <- data[!duplicated(data$`Component Name`), ]
  data <- data[order(data$`Retention time`), ]
  data[, "Component Name"] <- whto(a = LCdatabase[, c("Number", "Name")], b = data$`Component Name`)

  if (any(is.na(data[, "Component Name"]))) {
    stop("数据库不准确")
  }

  for (i in 1:dim(data)[1]) {
    label2 <- datalabel[datalabel$`Retention time` > data$`Retention time`[i], ]

    if (dim(label2)[1] == 0) {
      label2 <- datalabel[which.max(datalabel$`Retention time`), ]
    } else {
      label2 <- label2[which.min(label2$`Retention time`), ]
    }

    data[i, 3:dim(data)[2]] <- data[i, 3:dim(data)[2]] / label2[, 3:dim(data)[2]]
  }


  data[, "mode"] <- paste0("LC-", mode)
  data <- data[, c(1, 2, dim(data)[2], 4:dim(data)[2] - 1)]
  data <- data[, c(T, T, T, (colnames(data)[-1:-3] %in% name[, 1]))]
  colnames(data)[-1:-3] <- whto(name, colnames(data)[-1:-3])

  return(data)
}
