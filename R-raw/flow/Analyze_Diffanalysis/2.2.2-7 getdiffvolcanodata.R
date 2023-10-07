#!/opt/conda/bin/Rscript

#' 根据rds获取差异的热图数据，并保存
#'
#' @param rdspath rds路径
#' @param savepath 保存路径
#'
#' @export
getdiffvolcano_file <- function(rdspath,
                                savepath = NULL){
  
  # 数据提取
  data <- readdata(rdspath)
  
  volcanodata <- get_diffana_data(data = data,wb = NULL,
                                  datafrom = "diffdata",
                                  needlist = c("Metabolites","Accession"),
                                  sort = c("log2FoldChange","p-value","VIP"),
                                  decreasing = c(F,F,T))
  
  volcanodata <- volcanodata[,!(colnames(volcanodata) %in% row.names(data[["args"]][["singleclass"]])),drop = F]
  
  if(is.null(savepath)){
    
    return(volcanodata)
    
  }else{
    # 数据保存
    savetxt(data = volcanodata,
            filename = savepath)
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-r","--rdspath",help = "diff_filter中的rds文件",required = T)
  parser$add_argument("-s","--savepath", default = NULL,help="结果存储路径")
  
  args <- parser$parse_args()
  
  result <- do.call(what = getdiffvolcano_file, args = args)
  
}
