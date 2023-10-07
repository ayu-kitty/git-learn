#!/opt/conda/bin/Rscript

# FC计算
#' @export
fc_ana <- function(data,
                   class,
                   compare,
                   log = F){
  
  data1 <- data
  
  if(any(data < 0,na.rm = T)){
    log <- T
  }
  
  group <- unlist(strsplit(compare,split = "\\/"))
  
  data1[,"id"] <- row.names(data1)
  data1 <- reshape2::melt(data1,value.name = "value",variable.name="Sample",id.vars ="id")
  data1 <- merge(data1,class,by.x= "Sample",by.y=0)
  data1 <- reshape2::dcast(data = data1,formula = id~Group,value.var = "value",
                           fun.aggregate = mean,
                           na.rm = T)
  data1 <- data1[,c("id",group),drop = F]
  data1[is.na(data1)] <- NA
  
  if(length(group) == 2){
    if(log){
      for ( i in 1:dim(data1)[1]) {
        if(is.na(data1[i,group[1]]) & is.na(data1[i,group[2]])){
          data1[i,"log2FoldChange"] <- 0
        }else if(is.na(data1[i,group[1]])){
          data1[i,"log2FoldChange"] <- -15
        }else if(is.na(data1[i,group[2]])){
          data1[i,"log2FoldChange"] <- 15
        }else{
          data1[i,"log2FoldChange"] <- data1[i,group[1]]-data1[i,group[2]]
        }
      }
      # data1[,"log2FoldChange"] <- data1[,group[1]]-data1[,group[2]]
      data1[,"FoldChange"] <- 2^data1[,"log2FoldChange"]
      data1[!is.na(data1[,"log2FoldChange"]) & data1[,"log2FoldChange"] > 0,"Regulation"] <- "Up"
      data1[!is.na(data1[,"log2FoldChange"]) & data1[,"log2FoldChange"] < 0,"Regulation"] <- "Down"
    }else{
      for ( i in 1:dim(data1)[1]) {
        if(is.na(data1[i,group[1]]) & is.na(data1[i,group[2]])){
          data1[i,"FoldChange"] <- 1
        }else if(is.na(data1[i,group[1]])){
          data1[i,"FoldChange"] <- 1/32768
        }else if(is.na(data1[i,group[2]])){
          data1[i,"FoldChange"] <- 32768
        }else{
          if(data1[i,group[1]] == 0 & data1[i,group[2]] == 0){
            data1[i,"FoldChange"] <- 1
          }else if(data1[i,group[1]] == 0){
            data1[i,"FoldChange"] <- 1/32768
          }else if(data1[i,group[2]] == 0){
            data1[i,"FoldChange"] <- 32768
          }else{
            data1[,"FoldChange"] <- data1[,group[1]]/data1[,group[2]]
          }
        }
      }
      # data1[,"FoldChange"] <- data1[,group[1]]/data1[,group[2]]
      data1[,"log2FoldChange"] <- log2(data1[,"FoldChange"])
      data1[!is.na(data1[,"FoldChange"]) & data1[,"FoldChange"] > 1,"Regulation"] <- "Up"
      data1[!is.na(data1[,"FoldChange"]) & data1[,"FoldChange"] < 1,"Regulation"] <- "Down"
    }
  }
  
  for ( i in 1:length(group)) {
    colnames(data1)[colnames(data1) == group[i]] <- paste0("Average(",group[i],")")
  }
  diffdata <- data1 
  row.names(diffdata) <- diffdata[,"id"]
  diffdata <- diffdata[,colnames(diffdata) != "id"]
  # diffdata <- diffdata[order(as.numeric(row.names(diffdata))),]
  
  return(diffdata)
}
