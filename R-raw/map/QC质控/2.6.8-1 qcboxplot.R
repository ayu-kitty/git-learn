#!/opt/conda/bin/Rscript

#' @export
qc_boxplot <- function(data = c("数据矩阵.xlsx","数据矩阵"),
                       class1 = c("数据矩阵.xlsx","分组"),
                       log = T,
                       num = 60,
                       needgroup = NULL,
                       mapname = "Intensity Distribution",
                       savepath = "./",
                       imagetype=c("jpg","pdf"),
                       height = 6,
                       width = 5,
                       annotation_colors = c(
                         "green4", "blue3", "firebrick",
                                 "gold", "darkviolet", "darkorange",
                                 "skyblue3", "olivedrab3", "dodgerblue3",
                                 "aquamarine2", "deeppink3", "slateblue3",
                                 "brown2", "palegreen2", "chocolate2",
                                 "antiquewhite3", "steelblue1", "violetred1",
                                 "burlywood3", "pink1", "slategray2",
                                 "orangered1", "cyan3", "yellow4",
                                 "red", "plum", "greenyellow",
                                 "mediumpurple2", "tan1", "magenta"
                       ),
                       ...){
  suppressMessages(library("dplyr"))
  
  data <- readdata(filename = data)
  class1 <- readdata(filename = class1)
  
  ID<-data.frame(data$ID)
  data <- data[,colnames(data) %in% class1$sample]
  boxdata<-cbind(ID,data,row.names =T)
  
  if (log) {
    boxdata <- log10(boxdata)
  }
  
  if (dim(boxdata)[2] > num) {
    boxdata <- boxdata[, 1:num]
  }
  
  data1 <- boxdata[, 0]
  class <- NULL
  class1[2]%>%.[!duplicated(.),]%>%setdiff(.,"QC")%>%union("QC",.)->needgroup1
  
  for (k in 1:length(needgroup1)) {
    partdata <- boxdata[, colnames(boxdata) %in%
                          unlist(strsplit(class1$sample[class1$`Group-1` %in% needgroup1[k]],
                                          split = ","
                          )), drop = F]
    data1 <- cbind(data1, partdata)
    class<- c(class, rep(needgroup1[k], dim(partdata)[2]))
  }
  
  boxdata <- data1
  
  boxdata <- data.frame(
    metabolites = row.names(boxdata),
    boxdata,
    check.names = F, stringsAsFactors = F
  )
  boxdata <- reshape2::melt(boxdata, id = c("metabolites"))
  colnames(boxdata) <- c("Metabolites", "Sample", "Expression")
  boxdata[, "Group"] <- merge(boxdata[, "Sample", drop = F],
                              data.frame(
                                Sample = colnames(data1),
                                class = class, stringsAsFactors = F
                              ),
                              sort = F
  )[, 2]
  
  boxdata$Group <- factor(boxdata$Group,levels = unique(boxdata$Group))
  
  suppressMessages(library("ggplot2"))
  pp <- ggplot(boxdata, aes(x = Sample, y = Expression)) +
    theme_bw() +
    geom_boxplot(aes(fill = Group),outlier.size = ifelse((length(unique(class1$`Group-1`)) + 2) > 12, 0.5, (12 / length(unique(class1$`Group-1`)) ))) +
    theme(plot.title = ggplot2::element_text(hjust = 0.5),
          axis.text.y = ggplot2::element_text(size = 12),
          axis.title = ggplot2::element_text(size = 12),
          panel.grid.minor =element_blank(),
          panel.grid.major =element_blank(), 
          panel.background = element_blank(),
          # BUG 列名不存在Group-1怎么处理
          axis.text.x = element_text(angle = 45, hjust = 1, size = ifelse((length(unique(class1$`Group-1`)) + 2) > 12, 
                                                                          5, 
                                                                          (length(unique(class1$`Group-1`)))))) +
    labs(x = "sample", y = "log10(Intensity)",title = "Metabolites Intensity Distribution") +
    scale_fill_manual(values = annotation_colors)
  
  ggplotsave(plot = pp,
             savepath = saveptah,
             mapname = mapname,
             imagetype = imagetype,
             width = width,
             height = height,
             ...)
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_qc_boxplot <- map_autodraw$new(moudle = qc_boxplot)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "饼图矩阵文件",required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",
                      help = "xlsx中的sheet，全部分析请忽略")
  
  parser$add_argument("-cl","--class1",
                      default = c("数据矩阵.xlsx","分组"), 
                      nargs = "+",
                      help = "分组文件,如果是xlsx文件以`数据矩阵.xlsx 分组`形式传参",
                      required = T)
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf","html"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-l","--log",default = T, help = "是否log")
  parser$add_argument("-n","--num",default = 60, help = "样本数量")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  result <- do.call(what = map_qc_boxplot,args = args)
}

#' 根据文件进行饼图可视化
#' 
#' @export
map_qc_boxplot <- map_autodraw$new(moudle = qc_boxplot)$draw
