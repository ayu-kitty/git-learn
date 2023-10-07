#!/opt/conda/bin/Rscript

#' roc可视化
#'
#' @param data 数据
#' @param height 图片长度
#' @param width 图片宽度
#' @param legend.title 图例标题主题
#' @param plot.title 图标题主题
#' @param collist 颜色列表
#' @param col 颜色
#' @param aspect.ratio 长宽比
#' @param axis.text 轴文本主题
#' @param axis.title 轴标题主题
#' @param legend.text 图例文本主题
#' @param plot.margin 图空白边框
#' @param AUC 逻辑，是否显示auc值
#' @param cex auc值大小
#' @param hjust auc值对齐方向
#' @param nudge_x auc值x轴调整
#' @param nudge_y auc值y轴调整
#' @param right 逻辑，是否auc>0.5
#' @param legendcolor 图例标题
#' @param linesize 线粗细
#' @param linealpha 线透明度
#' @param label.angle 标签方向
#' @param ... 
#'
#' @export
auto_plotroc <- function(data,
                         collist = c(
                           "black", "green4", "blue3",
                           "firebrick", "gold", "darkviolet",
                           "darkorange", "skyblue3", "olivedrab3",
                           "deeppink3", "slateblue3", "brown2"
                         )[ifelse(length(unique(data$Feature)) == 1, 1, -1)],
                         col = collist[1:length(unique(data$Feature))],
                         aspect.ratio = 1,
                         axis.text = ggplot2::element_text(size = 10),
                         axis.title = ggplot2::element_text(size = 10),
                         legend.text = ggplot2::element_text(size = 6),
                         legend.title = ggplot2::element_text(size = 8),
                         plot.margin = ggplot2::unit(c(0.3, 0.3, 0.3, 0.3), "in"),
                         plot.title = ggplot2::element_text(hjust = 0.5),
                         AUC = T,
                         cex = 3,
                         hjust = 0,
                         vjust = 1,
                         nudge_x = 0,
                         nudge_y = 0,
                         right = T,
                         legendcolor = "Feature",
                         linesize = 0.6,
                         label.angle = 0,
                         linealpha = 0.8,
                         height = 5,
                         width = 5,
                         ...) {
  suppressMessages(library("ggplot2"))
  suppressMessages(library("plotROC"))
  suppressMessages(library("pROC"))

  data1 <- data

  pp <- ggplot(data1, aes(m = Expression,
                         d = Group,
                         color = Feature)) +
    geom_roc(labels = F, n.cuts = 0, size = linesize,linealpha = linealpha) +
    theme_bw()

  auc <- calc_auc(pp)[, c("Feature", "AUC")]
  names(auc)[1:2] <- c("Feature", "AUC")
  auc[, 1] <- unique(data1[, "Feature"])
  auc$adjAUC <- auc$AUC

  if(right){
    if(!all(auc$AUC >= 0.5)){
      pp <- ggplot(data1, aes(m = Expression,
                             d = Group,
                             color = Feature)) +
        geom_roc(mapping = aes(m = -Expression,
                               d = Group,
                               color = Feature),
                 data = data1[data1$Feature %in% auc[auc$AUC < 0.5, 1],],
                 inherit.aes = F,
                 linealpha = linealpha,
                 labels = F, n.cuts = 0, size = linesize)
      if(any(auc$AUC >= 0.5)){
        pp <- pp + geom_roc(mapping = aes(m = Expression,
                                       d = Group,
                                       color = Feature),
                         data = data1[data1$Feature %in% auc[auc$AUC >= 0.5, 1],],
                         inherit.aes = F,
                         linealpha = linealpha,
                         labels = F, n.cuts = 0, size = linesize)
      }

      pp <- pp+theme_bw()

      auc$adjAUC[auc$AUC < 0.5] <- 1-auc$adjAUC[auc$AUC < 0.5]
    }
  }

  if (AUC) {
    auctext <- round(auc$adjAUC, 3)
    if (length(unique(data$Feature)) == 1) {
      auctext <- paste0("AUC=", auctext)
      pp <- pp + annotate("text",
                       x = 0.5 + nudge_x,
                       y = 0.5 + nudge_y,
                       label = auctext,
                       angle = label.angle,
                       size = cex, 
                       hjust = hjust,
                       vjust = vjust)
    } else {
      auctext <- paste(unique(data$Feature), format(auc$adjAUC, digits = 5, nsmall = 3), sep = ":", collapse = "\n")
      auctext <- paste0("AUC","\n",auctext)
      pp <- pp + annotate("text",
                        x = 0.5 + nudge_x,
                        y = 0.5 + nudge_y,
                        label = auctext,
                        angle = label.angle,
                        size = cex, 
                        hjust = hjust,
                        vjust = vjust)
    }
  }

  if (length(unique(data$Feature)) == 1) {
      pp <- pp + facet_wrap(~Feature) +
        theme(legend.position = "none")
  }

  for (i in 1:dim(auc)[1]) {

    if(right & auc[i, "AUC"] < 0.5){
      a <- calculate_roc(M = -data1[data1$Feature == auc[i, 1], "Expression"],
                         D = data1[data1$Feature == auc[i, 1], "Group"],
                         ci = FALSE,
                         alpha = 0.05)
    }else{
      a <- calculate_roc(M = data1[data1$Feature == auc[i, 1], "Expression"],
                         D = data1[data1$Feature == auc[i, 1], "Group"],
                         ci = FALSE,
                         alpha = 0.05)
    }

    if (auc[i, "adjAUC"] < 0.5) {
      b <- a[which.min(1 - a$FPF + a$TPF), ]
    } else {
      b <- a[which.max(1 - a$FPF + a$TPF), ]
    }

    auc[i, "Specificity"] <- (1 - b[, 1])
    auc[i, "Sensitivity"] <- b[, 2]
    if(right & auc[i, "AUC"] < 0.5){
      auc[i, "cutoff"] <- -b[, 3]
    }else{
      auc[i, "cutoff"] <- b[, 3]
    }

    auc1 <- roc(data1[data1$Feature == auc[i, 1], "Group"], data1[data1$Feature == auc[i, 1], "Expression"], ci = T)
    if (auc[i, "adjAUC"] < 0.5) {
      auc[i, "95%ci"] <- paste0(
        format((1 - auc1[["ci"]][3]), digits = 3, nsmall = 3), "-",
        format((1 - auc1[["ci"]][1]), digits = 3, nsmall = 3)
      )
    } else {
      auc[i, "95%ci"] <- paste0(
        format(auc1[["ci"]][1], digits = 3, nsmall = 3), "-",
        format(auc1[["ci"]][3], digits = 3, nsmall = 3)
      )
    }
  }

  pp <- pp + xlab("1 - Specificity") +
    ylab("Sensitivity") +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    geom_hline(yintercept = c(0, 1),linetype = "dashed",alpha = 0.5) +
    geom_vline(xintercept = c(0, 1),linetype = "dashed",alpha = 0.5) +
    scale_color_manual(values = col, breaks = unique(data1[,"Feature"])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",alpha = 0.5) +
    labs(color = legendcolor) +
    theme(panel.grid = element_blank(),
          plot.title = plot.title,
          aspect.ratio = aspect.ratio,
          axis.text = axis.text,
          axis.title = axis.title,
          legend.text = legend.text,
          legend.title = legend.title,
          plot.margin = plot.margin)

  ggplotsave(plot = pp,
             height = height,
             width = width,
             ...)

  return(list(auc = auc,
              plot = pp))
}


#' roc可视化
#'
#' @param data 数据
#' @param mapname 保存文件名
#' @param number 图上绘制roc数据
#' @param saveroc 保存roc数据
#' @param ... 见[auto_plotroc()]
#'
#' @export
auto_roc <- function(data,
                     mapname = "ROC",
                     number = 1,
                     saveroc = T,
                     savepath = "./",
                     ...) {

  data1 <- data

  if (row.names(data1)[1] != "Group") {
    stop("请在数据第二行添加Group分组信息")
  }

  group <- table(as.vector(t(data1[1, ])))

  if (any(group < 3)) {
    savetxt(data = "有比较组分析未提供roc分析，可能由于分组样本小于3个或比较组大于两组",
            filename = paste0(savepath,"/说明.txt"))
    warning("有分组样本小于3个", immediate. = T)
    return("有分组样本小于3个")
  } else if (length(group) > 2) {
    savetxt(data = "有比较组分析未提供roc分析，可能由于分组样本小于3个或比较组大于两组",
            filename = paste0(savepath,"/说明.txt"))
    warning("比较组大于两组", immediate. = T)
    return("比较组大于两组")
  }

  if (is.null(number)) {
    data1 <- as.data.frame(t(data1))
    data1[, 2:dim(data1)[2]] <- apply(data1[, 2:dim(data1)[2], drop = F], 2, as.numeric)
    data1 <- reshape2::melt(data1, id = c("Group"))
    names(data1) <- c("Group", "Feature", "Expression")

    auc <- auto_plotroc(data = data1,
                        savepath = savepath,
                        mapname = mapname,
                        ...)

    auc$auc[, "ROC图名称"] <- mapname
  } else {
    auc1 <- NULL
    i <- 1
    len <- dim(data)[1]
    while ((i - 1) * number + 2 <= len) {
      rowdata <- c(1, (((i - 1) * number) + 2):(((i) * number) + 1))
      rowdata <- rowdata[rowdata %in% 1:dim(data1)[1]]
      data2 <- data1[rowdata, ]
      data2 <- as.data.frame(t(data2))
      data2[, 2:dim(data2)[2]] <- apply(data2[, 2:dim(data2)[2], drop = F], 2, as.numeric)
      data2 <- reshape2::melt(data2, id = c("Group"))
      names(data2) <- c("Group", "Feature", "Expression")

      auc <- auto_plotroc(data = data2,
                          mapname = paste0(mapname, "-", i),
                          savepath = savepath,
                          ...)
      
      # print(i)
      auc$auc[, "ROC图名称"] <- paste(mapname, "-", i, sep = "")
      i <- i + 1
      auc1 <- rbind(auc1, auc$auc)
    }

    auc[["auc"]] <- auc1
  }

  # print("auto_roc运行完成")
  
  if (saveroc) {
    savexlsx1(data = auc$auc, 
              filename = paste0(savepath,"/AUC.xlsx"), 
              sheet = mapname)
  } else {
    return(auc)
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_common_roc <- map_autodraw$new(moudle = auto_roc,row.names = 1)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "表达数据矩阵",
                      required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "分析结果", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf","html"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf","html"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 特有参数
  parser$add_argument("-nsr","--nosaveroc", default = T, action = "store_false", 
                      help="是否保存roc结果数据",dest = "saveroc")
  parser$add_argument("-a","--AUC", default = "T", help="是否显示AUC")
  parser$add_argument("-r","--right",default = "T",help = "AUC值是否自动换算成>0.5")
  parser$add_argument("-nu","--number", default = 1, type = "integer", 
                      help="多roc叠加，要与输入的数据量对应")
  parser$add_argument("-z","--zip",default = F, help = "是否压缩",action='store_true')
  
  args <- parser$parse_args()
  
  zip <- args$zip
  args$zip <- NULL
  
  rocresult <- do.call(what = map_common_roc, args = args)
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}

  if(zip){
    zip::zip(zipfile = "分析结果.zip",files = args$saveptah)
    unlink(list.files(pattern = "[^z][^i][^p]$",include.dirs = T,all.files = T), recursive = T)
  }
  
}

#' 根据文件进行roc可视化
#' 
#' @export
map_common_roc <- map_autodraw$new(moudle = auto_roc,row.names = 1)$draw
