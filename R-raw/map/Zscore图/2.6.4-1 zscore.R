#!/opt/conda/bin/Rscript

#' zscore图
#'
#' @param data 数据
#' @param height 图片长度
#' @param width 图片宽度
#' @param other 其他参数
#' @param aspect.ratio 长宽比
#' @param shape 点形状
#' @param linetype p筛选线类型
#' @param size 点大小
#' @param alpha 点透明度
#' @param annotation_colors 颜色列表
#' @param col 颜色
#' @param title 标题
#' @param x.title x轴标题
#' @param y.title y轴标题
#' @param axis.text 轴文本主题
#' @param axis.title 轴标题主题
#' @param legend.text  图例文本主题
#' @param legend.title  图例标题主题
#' @param xrange 
#' @param ... 
#'
#' @export
zscoremap <- function(data,
                      mapname = "Z-score",
                      savepath = "./",
                      width = 16,
                      height = 10,
                      size = 4,
                      shape = 1,
                      alpha = 1,
                      classfile = "classtype.xlsx",
                      col = stylefun_group(classfile = classfile,object = t(data)[,1],styletype = "fill"),
                      xrange = c(-3,3),
                      linetype = 3,
                      title = "",
                      x.title = "Z-score",
                      y.title = "Feature",
                      axis.text = ggplot2::element_text(size = ifelse((5 + 50 / dim(data)[1]) > 5.5, 14, (14 + 50 / dim(data)[1]))),
                      axis.title = ggplot2::element_text(size = 20),
                      legend.text = ggplot2::element_text(size = 20),
                      legend.title = ggplot2::element_text(size = 20),
                      aspect.ratio = 16 / 10,
                      ...) {
  suppressMessages(library("stringr"))
  suppressMessages(library("ggplot2"))
  
  options(warn = -1)
  
  data1 <- data
  if (row.names(data1)[1] != "Group") {stop("请在数据第二行添加Group分组信息")}
  
  if (dim(data1)[1] < 3 | dim(data1)[1] < 3) {
    savetxt(data = "数据少于3，未提供zscore图",
            filename = paste0(savepath,"/说明.txt"))
    warning("数据少于3，未提供zscore图", immediate. = T)
    return()
  }
  
  group <- table(as.vector(t(data1[1, ])))
  
  if (any(group < 3)) {
    savetxt(data = "有比较组分析未提供zscore图，可能由于分组样本小于3个",
            filename = paste0(savepath,"/说明.txt"))
    warning("有分组样本小于3个", immediate. = T)
    return()
  }
  
  for(k in 1:nrow(data)){
    row.names(data)[k] <- unlist(strsplit(split = ";\n",x = row.names(data)[k]))[1]
    row.names(data)[k] <- unlist(strsplit(split = "; ",x = row.names(data)[k]))[1]
    if(str_length(row.names(data)[k])>70){
      newname <- paste0(substring(row.names(data)[k],1,35),"...")
      i <- 2
      while (newname %in% row.names(data)) {
        newname <- paste0(substring(row.names(data)[k],1,35),"...-",i)
        i <- i+1
      }
      row.names(data)[k] <- newname
    }
  }
  
  data1 <- data
  data1 <- as.data.frame(t(data1))
  data1[, 2:dim(data1)[2]] <- apply(data1[, 2:dim(data1)[2]], 2, as.numeric)
  data1[, 2:dim(data1)[2]] <- scale(data1[, 2:dim(data1)[2]])
  
  data1 <- reshape2::melt(data1, id = c("Group"))
  names(data1) <- c("Group", "metabolites", "zscore")
  
  pp <- ggplot(data1, aes(x = zscore, y = metabolites)) +
    labs(x = x.title, y = y.title, title = title) +
    geom_point(aes(colour = Group), size = size, shape = shape, alpha = alpha,stroke = 1.5) +
    scale_colour_manual(values = col) +
    xlim(xrange) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust=0.5,size = 20),#title居中
          axis.line = element_line(colour = "black")) +
    theme(axis.text = axis.text,
          axis.title = axis.title,
          legend.text = legend.text,
          legend.title = legend.title,
          panel.grid.major.y = element_line(colour = "grey", linetype = linetype),
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(fill=NA,color="black", linewidth=1.5),
          aspect.ratio = aspect.ratio)
  
  ggplotsave(plot = pp,
             height = height,
             width = width,
             mapname= mapname,
             savepath = savepath,
             ...)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_common_zscore <- map_autodraw$new(moudle = zscoremap,row.names = 1)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "表达数据矩阵",
                      required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf","html"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf","html"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  parser$add_argument("-cf","--classfile",default = "classtype.xlsx", help = "模板")
  
  # 此图参数
  parser$add_argument("-ti","--title", default = "", type = "character", help="图片标题")
  parser$add_argument("-xt","--x.title", default = "Z-score", type = "character", help="图片x轴标签")
  parser$add_argument("-yt","--y.title", default = "Feature", type = "character", help="图片y轴标签")
  
  args <- parser$parse_args()

  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}

  heatmapargs <- do.call(what = map_common_zscore, args = args)
}

#' 根据文件进行zscore可视化
#' 
#' @export
map_common_zscore <- map_autodraw$new(moudle = zscoremap,row.names = 1)$draw
