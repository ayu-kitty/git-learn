#!/opt/conda/bin/Rscript

#' 根据rds获取差异的热图数据，并保存
#'
#' @param rdspath rds路径
#' @param savepath 保存路径
#'
#' @export
getdiffheatmap_file <- function(rdspath,
                                savepath = NULL,
                                addgroup = T,
                                top = NULL,
                                scale = F,
                                order = T){
  
  # 数据提取
  data <- readdata(rdspath)
  
  heatmapdata <- get_diffana_data(data = data,wb = NULL,
                                  datafrom = "filterdata",
                                  needlist = c("Metabolites","Accession"),
                                  sort = c("log2FoldChange","p-value","VIP"),
                                  decreasing = c(F,F,T))
  
  if(!is.null(top)){
    
    if(dim(heatmapdata)[1] > top){
      heatmapdata <- heatmapdata[1:top,,drop = F]
    }
    
  }
  
  infolist <- c("mz","Metabolites","ID","Accession")
  infodata <- heatmapdata
  
  for ( i in 1:length(infolist)) {
    if(infolist[i] %in% colnames(infodata)){
      infodata <- infodata[,infolist[i],drop = F]
    }
  }
  
  if(order){
    diffdata <- heatmapdata[,unlist(data[["args"]][["anaclass"]][rev(data[["args"]][["group2"]])]),drop=F]
  }else{
    diffdata <- heatmapdata[,row.names(data[["args"]][["singleclass"]]),drop=F]
  }
  
  if(scale){
    diffdata2 <- t(apply(diffdata,1,scale))
    colnames(diffdata2) <- colnames(diffdata)
    diffdata <- diffdata2
  }
  
  if(addgroup){
    group <- as.data.frame(t(data[["args"]][["singleclass"]]))
    group <- group[,colnames(diffdata),drop = F]
    diffdata <- rbind(group,diffdata)
    
    infodata1 <- infodata[1,,drop =F]
    infodata1 <- "Group"
    infodata <- rbind(infodata1,infodata)
  }

  heatmapdata <- cbind(infodata,diffdata)
  
  if(is.null(savepath)){
    
    row.names(heatmapdata) <- heatmapdata[,1]
    heatmapdata <- heatmapdata[,-1,drop = F]
    
    return(heatmapdata)
  }else{
    # 数据保存
    savetxt(data = heatmapdata,
            filename = savepath)
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-r","--rdspath",help = "diff_filter中的rds文件",required = T)
  parser$add_argument("-s","--savepath", default = NULL,help="结果存储路径")
  parser$add_argument("-a","--addgroup", default = "T", help="是否添加组别")
  parser$add_argument("-t","--top", default = NULL, type = "integer", help="取前几的差异")
  parser$add_argument("-sc","--scale", default = F, action = "store_true", help="取前几的差异")

  args <- parser$parse_args()
  
  args$addgroup <- as.logical(args$addgroup)
  
  result <- do.call(what = getdiffheatmap_file, args = args)
  
}
