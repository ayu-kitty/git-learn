#!/opt/conda/bin/Rscript

#' 绘制小提琴图 带误差线 带箱线图
#'
#' @param data 数据
#' @param mapname 保存文件名
#' @param imagetype 保存图片格式
#' @param height 图片长度
#' @param width 图片宽度
#' @param aspect.ratio 长宽比
#' @param annotation_colors 颜色列标
#' @param col 颜色
#' @param axis.text 轴文本主题
#' @param axis.title 轴标题主题
#' @param legend.text 图例文本主题
#' @param legend.title 图例标题主题
#' @param shrink 分面轴分布
#' @param scales 分面轴缩放
#' @param nrow 分面行数量
#' @param ncol 分面列数量
#' @param legend 图例标题
#' @param xlab x轴标题
#' @param ylab y轴标题
#' @param add 添加绘图内容，见[ggviolin()]
#' @param compare 逻辑，添加比较
#' @param method 比较方法
#' @param ... 
#'
#' @export
auto_violin <- function(data,
                        mapname = NA,
                        imagetype = NA,
                        classfile = "classtype.xlsx",
                        col = stylefun_group(classfile = classfile,object = data$Group,styletype = "fill"),
                        height = 5,
                        width = 5,
                        axis.text = ggplot2::element_text(size = 10),
                        axis.title = ggplot2::element_text(size = 10),
                        legend.text = ggplot2::element_text(size = 10),
                        legend.title = ggplot2::element_text(size = 10),
                        aspect.ratio = 1,
                        shrink = T,
                        scales = "free",
                        nrow = NULL,
                        ncol = NULL,
                        add = "boxplot",
                        compare = F,
                        method = "t.test",
                        var.equal = T,
                        legend = "none",
                        xlab = "Group",
                        ylab = "Expression",
                        add.params = list(width = 0.2),
                        color = "Group",
                        ...) {
  suppressMessages(library("ggplot2"))
  suppressMessages(library("ggpubr"))
  pp <- ggviolin(data,
                 x = "Group", y = "Expression", color = color, add = add, 
                 facet.by = "Metabolites", add.params = add.params,
                 ...) +
    scale_color_manual(values = col, guide = legend) +
    scale_fill_manual(values = col, guide = legend) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    facet_wrap(~Metabolites, scales = scales, shrink = shrink, nrow = nrow, ncol = ncol) +
    theme(panel.grid = element_blank(),
          axis.text = axis.text,
          axis.title = axis.title,
          legend.text = legend.text,
          legend.title = legend.title,
          aspect.ratio = aspect.ratio)
  
  if (compare) {
    group <- unique(data$Group)
    
    if (length(group) < 4) {
      group <- combn(c(group), 2)
      da <- list()
      for (i in 1:ncol(group)) da[[i]] <- group[, i]
      pp <- pp + stat_compare_means(method = method, hide.ns = T, label = "p.signif", 
                                    comparisons = da, label.y = max(data$Expression * 1.35),
                                    method.args = list(var.equal = var.equal)) +
        scale_y_continuous(expand = expansion(mult = .3))
    }
  }
  
  ggplotsave(plot = pp,
             mapname = mapname,
             imagetype = imagetype,
             height = height,
             width = width,
             ...)
}

#' auto_violin2
#'
#' 绘制小提琴图 带误差线 带箱线图
#'
#' @param data 数据
#' @param name 保存文件名
#' @param type 保存图片格式
#' @param height 图片长度
#' @param width 图片宽度
#' @param family 图片字体
#' @param dpi 图片分辨率
#' @param compression tiff格式压缩模式
#' @param text_theme 文本主题
#' @param title_theme 标题主题
#' @param other 其他参数
#' @param aspect.ratio 长宽比
#' @param annotation_colors 颜色列标
#' @param col 颜色
#' @param axis.text 轴文本主题
#' @param axis.title 轴标题主题
#' @param legend.text 图例文本主题
#' @param legend.title 图例标题主题
#' @param shrink 分面轴分布
#' @param scales 分面轴缩放
#' @param nrow 分面行数量
#' @param ncol 分面列数量
#' @param legend 图例标题
#' @param xlab x轴标题
#' @param ylab y轴标题
#' @param add 添加绘图内容，见[ggviolin()]
#' @param compare 逻辑，添加比较
#' @param method 比较方法
#'
#' @export
auto_violin2 <- function(data,
                         mapname = NA,
                         imagetype = NA,
                         annotation_colors = list(c("green4", "blue3", "firebrick",
                                                    "gold", "darkviolet", "darkorange",
                                                    "skyblue3", "olivedrab3", "dodgerblue3",
                                                    "aquamarine2", "deeppink3", "slateblue3",
                                                    "brown2", "palegreen2", "chocolate2",
                                                    "antiquewhite3", "steelblue1", "violetred1",
                                                    "burlywood3", "pink1", "slategray2",
                                                    "orangered1", "cyan3", "yellow4",
                                                    "red", "plum", "greenyellow",
                                                    "mediumpurple2", "tan1", "magenta")[1:length(unique(data$Group))]),
                         col = annotation_colors[[1]],
                         height = 5,
                         width = 5,
                         axis.text = ggplot2::element_text(size = 10),
                         axis.title = ggplot2::element_text(size = 10),
                         legend.text = ggplot2::element_text(size = 10),
                         legend.title = ggplot2::element_text(size = 10),
                         aspect.ratio = 1,
                         shrink = T,
                         scales = "free",
                         nrow = NULL,
                         ncol = NULL,
                         compare = F,
                         method = "t.test",
                         var.equal = T,
                         legend = "none",
                         xlab = "Group",
                         ylab = "Expression",
                         color = "Group",
                         isscale = T,
                         ...) {
  
  if(isscale){
    data$Expression <- scale(data$Expression)
  }
  
  suppressMessages(library("ggplot2"))
  suppressMessages(library("ggpubr"))
  
  pp <- ggviolin(data,
                 x = "Group", y = "Expression", color = color,fill = color,
                 facet.by = "Metabolites",
                 ...) +
    scale_color_manual(values = col, guide = legend) +
    scale_fill_manual(values = col, guide = legend) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    facet_wrap(~Metabolites, scales = scales, shrink = shrink, nrow = nrow, ncol = ncol) +
    theme(panel.grid = element_blank(),
          axis.text = axis.text,
          axis.title = axis.title,
          legend.text = legend.text,
          legend.title = legend.title,
          aspect.ratio = aspect.ratio)
  
  if(isscale){
    pp <- pp + ylab("Normalized Expression")
  }
  
  pp <- pp+geom_boxplot(width=0.15,color="blue",
                        outlier.color = NA,size=0.2,mapping = aes(fill = Group))
  
  ggplotsave(plot = pp,
             mapname = mapname,
             imagetype = imagetype,
             height = height,
             width = width,
             ...)
}


#' auto_boxchart
#'
#' 箱线图类图像可视化
#'
#' @param data 数据
#' @param name 保存文件名
#' @param type 保存图片格式
#' @param dealname 是否处理行名
#' @param number 图上绘制数据
#' @param boxmoudle 绘图函数
#' @param ... boxmoudle函数参数
#'
#' @export
auto_boxchart <- function(data,
                          mapname = NULL,
                          imagetype = c("jpg", "pdf"),
                          number = 1,
                          boxmoudle = auto_violin,
                          savepath = "./",
                          ...) {
  
  suppressMessages(library("stringr"))
  
  group <- table(as.vector(t(data[1, ])))
  
  if (any(group < 3)) {
    savetxt(data = "有比较组分析未提供箱线图，可能由于分组样本小于3个",
            filename = paste0(savepath,"/说明.txt"))
    
    warning("有分组样本小于3个", immediate. = T)
    # return("有分组样本小于3个")
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
  
  if (is.null(number)) {
    data1 <- as.data.frame(t(data1))
    data1[, 2:dim(data1)[2]] <- apply(data1[, 2:dim(data1)[2]], 2, as.numeric)
    # data1[,2:dim(data1)[2]]<-scale(data1[,2:dim(data1)[2]])
    data1 <- reshape2::melt(data1, id = c("Group"))
    names(data1) <- c("Group", "Metabolites", "Expression")
    data1$Group <- factor(x = data1$Group, levels = unique(data1$Group))
    
    for (j in 1:length(type)) {
      returndata <- boxmoudle(
        data = data1,
        name = ifelse(is.na(type[j]), type[j], paste(name1, ".", type[j], sep = "")),
        type = type[j],
        savepath = savepath,
        ...
      )
    }
  } else {
    i <- 1
    len <- dim(data)[1]
    while ((i - 1) * number + 2 <= len) {
      data2 <- data1[c(1, (((i - 1) * number) + 2):(((i) * number) + 1)), ]
      data2 <- as.data.frame(t(data2))
      data2[, 2:dim(data2)[2]] <- apply(data2[, 2:dim(data2)[2], drop = F], 2, as.numeric)
      # data2[,2:dim(data2)[2]]<-scale(data2[,2:dim(data2)[2],drop=F])
      data2 <- reshape2::melt(data2, id = c("Group"))
      names(data2) <- c("Group", "Metabolites", "Expression")
      data2$Group <- factor(x = data2$Group, levels = unique(data2$Group))
      
      returndata <- boxmoudle(data = data2,
                              mapname = paste0(mapname,"-",i),
                              imagetype = imagetype,
                              savepath = savepath,
                              ...)
      
      print(i)
      i <- i + 1
    }
  }
  
  print("auto_boxchart运行完成")
  return()
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_common_boxplot <- map_autodraw$new(moudle = auto_boxchart,row.names = 1)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "矩阵文件",
                      required = T)
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
  parser$add_argument("-co","--compare", default= F, action = "store_true", help = "是否显示显著性")
  parser$add_argument("-as","--aspect.ratio", default= 1, type= "double", help = "图片横纵比")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  heatmapargs <- do.call(what = map_common_boxplot, args = args)
}

#' @export
map_common_boxplot <- map_autodraw$new(moudle = auto_boxchart,row.names = 1)$draw

