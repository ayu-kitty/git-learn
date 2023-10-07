#!/opt/conda/bin/Rscript

#' 根据变量进行多元统计score图可视化
#' 
#' @param data 数据
#' @param mapname 图片名称
#' @param height 图片长度
#' @param width 图片宽度
#' @param x x轴标签
#' @param y y轴标签
#' @param title 标题
#' @param showname 逻辑，是否显示样本名
#' @param max.overlaps 名称重叠度
#' @param fill 点填充颜色
#' @param colour 点边框颜色
#' @param shape 点形状
#' @param other 其他参数
#' @param ... 见`ggplotsave`
#'
#' @export
mulstatisticsmap2 <- function(data,
                              mapname = "statistics",
                              x = colnames(data)[2],
                              y = colnames(data)[3],
                              mode = "",
                              title = mode,
                              classfile = "classtype.xlsx",
                              point_fill = stylefun_group(classfile = classfile,object = data,styletype = "fill"),
                              point_colour = stylefun_group(classfile = classfile,object = data,styletype = "colour"),
                              point_shape = stylefun_group(classfile = classfile,object = data,styletype = "shape"),
                              showname = F,
                              point_size=4,
                              fill = T,
                              color = F,
                              max.overlaps = 200,
                              other = NULL,
                              width = 8,
                              height = 5,
                              imagetype = c("jpg","pdf","html"),
                              ...) {
  suppressMessages(library("ggplot2"))
  options(warn=-1)
  
  score <- data
  names(score)[2] <- c("p1")
  names(score)[3] <- c("p2")
  
  score$Group <- factor(score$Group, levels = unique(score$Group))
  
  centroids <- aggregate(cbind(p1,p2)~ Group,score,mean)
  row.names(centroids) <- centroids$Group
  conf.rgn <- do.call(rbind,lapply(unique(score$Group),function(t){
    if(sum(score$Group==t) == 1){
      NULL
    }else{
      data.frame(Group=as.character(t),
                 ellipse::ellipse(cov(score[score$Group==t,2:3]),
                                  centre=as.matrix(centroids[t,2:3]),level=0.95),
                 stringsAsFactors=FALSE)
    }
  }))
  pp <- ggplot(data = score, mapping = aes(x = p1, 
                                           y = p2))
  if(fill == T & color == F){
    pp <- pp + geom_polygon(data = conf.rgn,
                 alpha = 0.3,
                 mapping = aes(fill=Group),
                 show.legend = F)}else if(color == T & fill == F){
    pp <- pp + geom_polygon(data = conf.rgn,
                 alpha = 0.3,
                 mapping =  aes(color=Group),fill=NA,
                 show.legend = F)}else if(color == T & fill == T){
    pp <- pp + geom_polygon(data = conf.rgn,
                 alpha = 0.3,
                 mapping =  aes(color=Group,fill=Group),
                 show.legend = F)}
  pp <- pp+labs(title = title, x = x, y = y)+
    geom_hline(yintercept = 0,linetype = 2,col = "grey")+
    geom_vline(xintercept = 0,linetype = 2,col = "grey")+
    geom_point(size = point_size,
               mapping = aes(shape = Group,
                             fill = Group,
                             colour = Group,
                             text = paste0("Sample:",Sample, "<br>Group:", Group)))+
    # stat_ellipse(map = aes(fill = Group),
    #              geom = "polygon",
    #              type = "norm",
    #              alpha = 0.3,
    #              segments = 200, 
    #              show.legend = F)+
    scale_shape_manual(values = point_shape)+
    scale_fill_manual(values = point_fill)+
    scale_colour_manual(values = point_colour)+
    guides(fill="none")+
    guides(shape="none")+
    theme_bw()+
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 15),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black",linewidth = 1),
          # axis.line = element_line(colour = "black"),
          aspect.ratio = 3/4)+
    expand_limits(y = c(-max(abs(fivenum(score$p2)))*1.3,
                        max(abs(fivenum(score$p2)))*1.3),
                  x = c(-max(abs(fivenum(score$p1)))*1.3,
                        max(abs(fivenum(score$p1)))*1.3))

  if (showname) {
    pp <- pp + ggrepel::geom_text_repel(mapping = aes(label = Sample),
                                        max.overlaps = max.overlaps,
                                        na.rm = T)
  }
  
  args <- ggplotsave(plot = pp,
                     mapname = mapname,
                     other = other,
                     width = width,
                     height = height,
                     imagetype = imagetype,
                     ...)
  return(args)
}

#' @export
mulstatisticsmap2_obj <- function(data,...){
  if("mulstatistics" %in% class(data)){
    data <- getscore_file(rdspath = data)
  }
  result <- mulstatisticsmap2(data = data,...)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_mulstatistics_scoremap2 <- map_autodraw$new(mulstatisticsmap2_obj)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "score矩阵文件",required = T)
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
  
  # 此图参数
  parser$add_argument("-c","--classfile",default = "classtype.xlsx", help = "模板")
  parser$add_argument("-fi","--point_fill",default = "", help = "填充颜色",nargs="+")
  parser$add_argument("-co","--point_colour",default = "", help = "点外圈颜色",nargs="+")
  parser$add_argument("-ps","--point_shape",default = "", help = "点形状",nargs="+")
  parser$add_argument("-pi","--point_size",default = 4, type= "double",help = "点大小")
  parser$add_argument("-sn","--showname",default = F, type="logical",help = "是否显示样本名")
  parser$add_argument("-mo","--max_overlaps",default = 200, type= "double",help = "样本名重叠度",dest="max.overlaps")
  parser$add_argument("-md","--mode",default = "")
  parser$add_argument("-fl","--fill",default = T, type= "logical",help = "是否填充置信区间颜色")
  parser$add_argument("-col","--color",default = F, type= "logical",help = "置信区间边框颜色")
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}

  if("" %in% args$point_fill){args$point_fill <- NULL}
  if("" %in% args$point_colour){args$point_colour <- NULL}
  if("" %in% args$point_shape){args$point_shape <- NULL}
  
  result <- do.call(what = map_mulstatistics_scoremap2,args = args) 
}

#' 根据文件进行多元统计score图可视化
#' 
#' @export
map_mulstatistics_scoremap2 <- map_autodraw$new(mulstatisticsmap2_obj)$draw
