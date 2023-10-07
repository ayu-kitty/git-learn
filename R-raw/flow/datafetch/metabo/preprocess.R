#' preprocess
#'
#' 数据信息与数据分开存储
#'
#' @param data obj
#' @param order 样本名排序
#' @param del 逻辑，是否删除
#' @param mean 逻辑，是否均值
#'
#' @export
preprocess <- function(data,
                       order = T,
                       del = T,
                       mean = T) {
  print("数据预处理开始")

  data1 <- data
  
  if(!dir.exists("background")){
    movebackground(db = "background",
                   dbfrom = "/data/hstore4/database/kegg/kegg/",
                   keggmapfrom = "/data/hstore4/database/kegg/kegg/",
                   org = data[["info"]][["basic"]][["species"]],
                   inputfile = data$data$rawdata$data)
  }
  raw_data <- data$data$rawdata$data
  if(file.exists("./background/gene_reactome.backgroud.xls")){
    reacinfo <- readdata(filename="./background/gene_reactome.backgroud.xls",col.names=c("ChEBI","Reactome_id","Reactome_Description"),header=F)
    reacinfo$ChEBI<-gsub("ChEBI-","",reacinfo$ChEBI)
    new_data <- merge(reacinfo,raw_data,by="ChEBI",sort=F,all.y=T)
    reacnum <- which(colnames(raw_data)=="ChEBI")
    raw_data <- new_data[,c(4:(reacnum+3),1:3,(reacnum+4):ncol(new_data))]
  }
  kegginfo <- readdata(filename="./background/gene_kegg.backgroud.xls",col.names=c("KEGG","ID Annotation","Annotation"),header=F)
  new_data <- merge(kegginfo,raw_data,by="KEGG",sort=F,all.y=T)
  keggnum <- which(colnames(raw_data)=="KEGG")
  raw_data <- new_data[,c(4:(keggnum+3),1:3,(keggnum+4):ncol(new_data))]
  
  if("Score" %in% colnames(raw_data)){
    data$data$rawdata$data<-raw_data[order(-raw_data$Score),]
  }else{
    data$data$rawdata$data<-raw_data[order(raw_data$ID),]
  }
  
  information <- data$data$rawdata$data[, !colnames(data$data$rawdata$data) %in% data$info$sample$samplename, drop = F]
  row.names(information) <- information[, data$info$basic$索引$ID]

  data2 <- data$data$rawdata$data[, colnames(data$data$rawdata$data) %in% data$info$sample$samplename, drop = F]
  row.names(data2) <- information[, data$info$basic$索引$ID]

  if (del) {
    data2 <- data2[, !(whto(
      data.frame(
        samplename = data$info$sample$samplename,
        deal = data$info$sample$deal,
        stringsAsFactors = F
      ),
      colnames(data2)
    ) %in% "删除"), drop = F]
  }

  if (order) {
    data3 <- data2[, !(whto(
      data.frame(
        samplename = data$info$sample$samplename,
        class = data$info$sample$class[[1]], stringsAsFactors = F
      ),
      colnames(data2)
    ) %in% "QC"), drop = F]
    data3 <- data3[, order(whto(
      data.frame(
        samplename = data$info$sample$samplename,
        class = data$info$sample$class[[1]], stringsAsFactors = F
      ),
      colnames(data3)
    )), drop = F]
    data4 <- data2[, (whto(
      data.frame(
        samplename = data$info$sample$samplename,
        class = data$info$sample$class[[1]], stringsAsFactors = F
      ),
      colnames(data2)
    ) %in% "QC"), drop = F]
    data2 <- cbind(data4, data3)
  }

  if (mean) {
    i <- 1
    while (i <= dim(data2)[2]) {
      if (!is.na(data$info$sample$deal[grepl(paste0("^", colnames(data2)[i], "$"), data$info$sample$samplename)])) {
        data2[, i] <- apply(data2[, whto(
          data.frame(
            samplename = data$info$sample$samplename,
            class = data$info$sample$class[[1]],
            stringsAsFactors = F
          ),
          colnames(data2)
        ) %in% whto(
          data.frame(
            samplename = data$info$sample$samplename,
            class = data$info$sample$class[[1]],
            stringsAsFactors = F
          ),
          colnames(data2)
        )[i], drop = F],
        MARGIN = 1, FUN = mean
        )
      }
      i <- i + 1
    }
  }

  data1$data$predata$information <- information
  data1$data$predata$data <- data2
  data1$data$data$data <- data2

  print("数据预处理完毕")
  return(data1)
}
