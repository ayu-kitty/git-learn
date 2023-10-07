#!/opt/conda/bin/Rscript

chart.Correlation2 <- function(R,
                               histogram = TRUE,
                               method = c("pearson", "kendall", "spearman"),
                               ...) {
  x <- checkData(R, method = "matrix")
  if (missing(method)) {
    method <- method[1]
  }
  cormeth <- method 
  panel.cor <- function(x, y,
                        digits = 2, prefix = "", use = "pairwise.complete.obs",
                        method = cormeth, cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) {
      cex <- 0.6 / strwidth(txt)
    }
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    Signif <- symnum(test$p.value,
                     corr = FALSE, na = FALSE,
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " ")
    )
    text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3) / 1.3)
    text(0.5, 0.8, Signif, cex = cex, col = 2)
  }
  
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  dotargs <- list(...)
  dotargs$method <- NULL
  rm(method)
  hist.panel <- function(x, ... = NULL) {
    par(new = TRUE)
    hist(x,
         col = "light gray", probability = TRUE, axes = FALSE,
         main = "", breaks = "FD"
    )
    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
    rug(x)
  }
  if (histogram) {
    pairs(x,
          gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor,
          diag.panel = hist.panel
    )
  } else {
    pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor)
  }
}


#' @export
chartCorrelation <- function(data,
                             savepath = "./",
                             mapname = "Correlation",
                             imagetype = c("jpg","pdf"),
                             dpi = 300,
                             width = 5,
                             height = 5,
                             family = "sans",
                             units = "in",
                             ...) {
  plotfile(savepath = savepath,
           mapname = mapname,
           imagetype = imagetype,
           height = height,
           width = width,
           dpi = dpi,
           family = family,
           units = units)
  
  suppressMessages(library("PerformanceAnalytics"))

  chart.Correlation2(data, ...)
  
  plotsave()
}


#' @export
auto_chart.Correlation <- function(data,
                                   mapname = "Correlation",
                                   savepath = "./",
                                   ...) {
  data1 <- data
  
  if (dim(data1)[1] == 0) {
    
    savetxt(data = "数据为空,不进行相关性分析",
            filename = paste0(savepath,"/说明.txt"))
    
    return("无数据进行相关性分析")
  }
  

  chartCorrelation(data = data1,
                   mapname = mapname,
                   savepath = savepath,
                   ...)
  
  # print("auto_chart.Correlation运行完成")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_common_corrplot3 <- map_autodraw$new(moudle = auto_chart.Correlation,row.names = 1)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "x表达矩阵文件",required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  parser$add_argument("-fy","--filenamey",default = NULL,nargs="*",
                      help = "y表达矩阵文件,如果是xlsx文件以`数据矩阵.xlsx 数据矩阵`形式传参")
  
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
  parser$add_argument("-cm","--method",default = "pearson", help = "相关性计算方法,包括pearson,spearman,kendall")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  corrdata <- do.call(what = map_common_corrplot3,args = args)
  
}

#' 根据文件进行Correlation可视化
#' 
#' @export
map_common_corrplot3 <- map_autodraw$new(moudle = auto_chart.Correlation,row.names = 1)$draw
