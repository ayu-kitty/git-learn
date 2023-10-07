#!/opt/conda/bin/Rscript

#' 根据变量进行多元统计Loading图可视化
#' 
#' @param data 数据
#' @param mapname 图片名称
#' @param height 图片长度
#' @param width 图片宽度
#' @param x x轴标签
#' @param y y轴标签
#' @param title 标题
#' @param other 其他参数
#' @param ... 见`ggplotsave`
#'
#' @export
loadingmap2 <- function(data,
                        mapname = "Loading",
                        x = colnames(data)[2],
                        y = colnames(data)[3],
                        title = "Loading",
                        other = NULL,
                        width = 8,
                        height = 5,
                        imagetype = c("jpg","pdf","html"),
                        point_fill = c("#c4c4c4","#c0d09d"),
                        point_colour =c("#6b6b6b","#446a37"),
                        point_shape = 21,
                        point_size = c(1,4),
                        Featurelist = "ID",
                        ...) {
  
  suppressMessages(library("ggplot2"))
  options(warn=-1)
  
  score <- data
  names(score)[2] <- c("p1")
  names(score)[3] <- c("p2")
  score[,"Feature"] <- apply(score[,Featurelist,drop=F],1,paste,collapse =";")
  score[,"colour"] <- "VIP < 1"
  score[score[,"VIP"] >= 1,"colour"] <- "VIP >= 1"
  
  pp <- ggplot(data = score) +
    geom_hline(yintercept = 0,linetype = 2,col = "grey")+
    geom_vline(xintercept = 0,linetype = 2,col = "grey")+
    geom_point(mapping = aes(x = p1, y = p2,
                             text = paste0("ID:", Feature),
                             size = VIP,
                             fill = colour,
                             colour = colour),
               shape = point_shape,
               show.legend = T)+
    scale_fill_manual(values = point_fill)+
    scale_size_continuous(range = point_size)+
    scale_colour_manual(values = point_colour)+
    labs(title = title, x = x, y = y)+
    theme_bw()+
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 15),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black",linewidth = 1),
          # axis.line = element_line(colour = "black"),
          aspect.ratio = 3/4)+
    expand_limits(y = c(-max(abs(fivenum(score$p2))) - 0.05,
                        max(abs(fivenum(score$p2))) + 0.05),
                  x = c(-max(abs(fivenum(score$p1))) - 0.05,
                        max(abs(fivenum(score$p1))) + 0.05))
  
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
loadingmap2_obj <- function(data,infofile = NULL,...){
  if("mulstatistics" %in% class(data)){
    data <- getloading_file(rdspath = data,infofile = infofile)
  }
  result <- loadingmap2(data = data,...)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_mulstatistics_loadingmap2 <- map_autodraw$new(loadingmap2_obj)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                        help = "loading矩阵文件",required = T)
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
  parser$add_argument("-fl","--Featurelist",default = "ID", help = "特征名的列名",nargs="+")
  parser$add_argument("-pf","--point_fill",default = c("#c4c4c4","#c0d09d"), help = "填充颜色",nargs=2)
  parser$add_argument("-pc","--point_colour",default = c("#6b6b6b","#446a37"), help = "点外圈颜色",nargs=2)
  parser$add_argument("-ps","--point_shape",default = 21, type= "integer", help = "点形状")
  parser$add_argument("-pz","--point_size",default = c(1,4),type= "integer",help = "点大小",nargs=2)
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}

  result <- do.call(what = map_mulstatistics_loadingmap2,args = args) 
}


#' 根据文件进行多元统计Loading图可视化
#' 
#' @export
map_mulstatistics_loadingmap2 <- map_autodraw$new(loadingmap2_obj)$draw
