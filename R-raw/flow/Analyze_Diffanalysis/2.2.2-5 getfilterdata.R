#!/opt/conda/bin/Rscript

#' @export
all_get_diff_filter_data <- function(savefilename = "差异代谢物.xlsx",
                                     datafrom = "filterdata",
                                     needlist = c("Metabolites","Accession"),
                                     ...){
  
  result <- all_get_diff_ana_data(savefilename = savefilename,
                                  datafrom =  datafrom,
                                  needlist = needlist,
                                  ...)
  
  return(result)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-rp","--rdspath",default = "oecloud/diff_filter/", help = "rds数据路径")
  parser$add_argument("-r","--rdsname",
                      default = NULL, 
                      help = "中间过程数据位置",
                      nargs = "+")
  parser$add_argument("-s","--savefilename",default = "差异代谢物.xlsx", help = "保存结果")
  
  args <- parser$parse_args()
  
  if(is.null(args$rdsname)){
    args$rdsname <- NULL
  }
  
  flowargs <- do.call(what = all_get_diff_filter_data,args = args)
  
}
