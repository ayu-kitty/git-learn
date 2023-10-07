#!/opt/conda/bin/Rscript

#' @export
QC_corrplotfun <- function(data,
                           log = T,
                           qc = T,
                           classfile = "classtype.xlsx",
                           width = 10,
                           height = 10,
                           mapname = "QCcor",
                           ...){
  
  suppressWarnings(library("ggplot2"))
  
  ggscatter <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    df <- data.frame(x = x, y = y)
    sp1 <- ggplot(df, aes(x=x, y=y)) +
      geom_point(col = '#B09C8566',alpha=.4) +
      geom_abline(intercept = 0, slope = 1, col = '#CD3333')
    return(sp1)
  }
  
  ggdehist <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    df <- data.frame(x = x)
    dh1 <- ggplot(df, aes(x=x)) +
      geom_histogram(aes(y=after_stat(density)), bins = 50, fill = 'steelblue', color='lightgray', alpha=.4) +
      geom_density(aes(y=after_stat(density))) + 
      theme_minimal()
    return(dh1)
  }
  
  if(qc){
    classdata <- readdata(paste0(dirname(classfile), "/classfile.yaml"))
    data <- data[,colnames(data) %in% classdata[["QC"]],drop = F]
  }
  
  num <- dim(data)[2]
  width <- width/10*num
  height <- height/10*num
  
  if(log){
    data2 <- log2(data+1)
  }else{
    data2 <- data
  }
  
  pp <- GGally::ggpairs(data2,
                        lower = list(continuous = GGally::wrap(ggscatter,size= 4)),
                        diag = list(continuous = GGally::wrap(ggdehist))) + 
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill=NA),
          axis.text =  element_text(color='black'))
  
  ggplotsave(plot = pp,
             width = width,
             height = height,
             mapname = mapname,
             ...)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_common_corrplot4 <- map_autodraw$new(moudle = QC_corrplotfun,row.names = 1)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "x表达矩阵文件",required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-cf","--classfile",default = "classtype.xlsx", help = "模板")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  corrdata <- do.call(what = map_common_corrplot4,args = args)
  
}

#' @export
map_common_corrplot4 <- map_autodraw$new(moudle = QC_corrplotfun,row.names = 1)$draw
