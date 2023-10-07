#!/opt/conda/bin/Rscript

lastC <- function(x) {
  y <- sub(" +$", "", x)
  p1 <- nchar(y)
  cc <- substr(y, p1, p1)
  return(cc)
}

orderPvalue <- function(treatment, means, alpha, pvalue, console) {
  n <- length(means)
  z <- data.frame(treatment, means)
  letras <- c(letters[1:26], LETTERS[1:26], 1:9, c(
    ".", "+",
    "-", "*", "/", "#", "$", "%", "&", "^", "[", "]", ":",
    "@", ";", "_", "?", "!", "=", "#", rep(" ", 2000)
  ))
  w <- z[order(z[, 2], decreasing = TRUE), ]
  M <- rep("", n)
  k <- 1
  k1 <- 0
  j <- 1
  i <- 1
  cambio <- n
  cambio1 <- 0
  chequeo <- 0
  M[1] <- letras[k]
  q <- as.numeric(rownames(w))
  while (j < n) {
    chequeo <- chequeo + 1
    if (chequeo > n) {
      break
    }
    for (i in j:n) {
      s <- pvalue[q[i], q[j]] > alpha
      if (s) {
        if (lastC(M[i]) != letras[k]) {
          M[i] <- paste(M[i], letras[k], sep = "")
        }
      } else {
        k <- k + 1
        cambio <- i
        cambio1 <- 0
        ja <- j
        for (jj in cambio:n) {
          M[jj] <- paste(M[jj], "",
                         sep = ""
          )
        }
        M[cambio] <- paste(M[cambio], letras[k], sep = "")
        for (v in ja:cambio) {
          if (pvalue[q[v], q[cambio]] <= alpha) {
            j <- j + 1
            cambio1 <- 1
          } else {
            break
          }
        }
        break
      }
    }
    if (cambio1 == 0) {
      j <- j + 1
    }
  }
  w <- data.frame(w, stat = M)
  trt <- as.character(w$treatment)
  means <- as.numeric(w$means)
  output <- data.frame(means, groups = M)
  rownames(output) <- trt
  if (k > 81) {
    cat("\n", k, "groups are estimated.The number of groups exceeded the maximum of 81 labels. change to group=FALSE.\n")
  }
  invisible(output)
}

#' 绘制boxplot 带显著性标记 可多组
#'
#' @param data 表格数据数据
#' @param mapname 保存文件名
#' @param imagetype 保存图片格式
#' @param height 图片长度
#' @param width 图片宽度
#' @param dpi 分辨率
#' @param family 字体
#' @param classfile 分组颜色表格
#' @param col 颜色
#' @param x.title x轴标题
#' @param y.title y轴标题
#' @param legend.text.size 图例文本字体大小
#' @param legend.title.size 图例标题字体大小
#' @param savepath 保存路径
#' @param grouporder 分组排序
#' @param x.angle x轴标题角度
#' @param axis.text.size 轴文本字体大小
#' @param axis.title.size 轴标题字体大小
#' @param title.size 标题字体大小
#' @param ... 
#'
#' @export
boxplot_new <- function(data="violin_top50-CG24-vs-NU0.xls",
                        mapname="Boxplot",
                        dpi=300,
                        x.title="",
                        y.title="Intensity",
                        height=NULL,
                        width=NULL,
                        family="sans",
                        imagetype=c("png","pdf"),
                        grouporder="NA",
                        classfile = "classtype.xlsx",
                        col = stylefun_group(classfile = classfile,styletype = "fill"),
                        savepath="./",
                        x.angle=45,
                        axis.text.size=11,
                        axis.title.size=12,
                        title.size=12,
                        legend.title.size=11,
                        legend.text.size=10){
  
  pacman::p_load(ggpubr,ggplot2,stringr)
  
  data1 <- readdata(data,row.names = 1)
  data2 <- as.data.frame(t(data1))
  out <- data2
  alpha_index <- colnames(out)
  
  if (grouporder != "NA") {
    custom_level <- unlist(strsplit(grouporder, split = ","))
  } else {
    custom_level <- unique(out$Group)
  }
  
  group_num <- length(unique(as.character((out$Group))))
  
  if(is.null(width)){
    if (group_num < 3) {
      width <- 1.8 + 1.0 * group_num
    } else {
      width <- 1.8 + 0.6 * group_num
    }
  }
  if(is.null(height)){
    height <- 3.6 + 0.2 * group_num
  }
  
  for (i in 2:(length(alpha_index))) {
    j <- i-1
    
    if(str_length(alpha_index[i])>40){
      midlength <- round(str_length(alpha_index[i])/2)
      alpha_index[i] <- paste0(substr(alpha_index[i],1,midlength),"\n",
                               substr(alpha_index[i],midlength+1,length(alpha_index)))
    }
    
    plotdata <- data.frame(Group = out$Group, abundance = out[, i])
    plotdata$abundance <- as.numeric(plotdata$abundance)
    plotdata$Group <- factor(plotdata$Group, levels = custom_level)
    res <- pairwise.t.test(plotdata$abundance, plotdata$Group,p.adjust.method = "none")
    
    # if(length(unique(plotdata$Group))>2){
    #   sink(file=paste0(savepath, mapname,"-", j, ".xls"))
    #   print(res)
    #   sink()
    # }
    
    pValue <- res$p.value
    pValue[is.nan(pValue)] <- 1
    pValue <- pValue[lower.tri(pValue, diag=TRUE)]
    
    my_annotations <- as.character(
      symnum(pValue,cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("***", "**", "*", ""))
    )
    my_annotations1 <- my_annotations[my_annotations != ""]
    
    wtq <- levels(plotdata$Group)
    lis <- combn(wtq, 2)
    my_comparisons <- tapply(lis, rep(1:ncol(lis), each = nrow(lis)), function(i) i)
    my_comparisons1 <- my_comparisons[my_annotations != ""]
    
    my_levels <- custom_level
    
    p <- ggplot(data = plotdata, aes(x = Group, y = abundance), colour = Group)+
      geom_boxplot(aes(fill=Group),alpha = 0.5,size = 0.6,width = 0.7,outlier.alpha = 0) +
      geom_jitter(aes(color = Group,size=4),alpha = 0.3,size = 1) +
      scale_color_manual(limits = my_levels,values = col) +
      scale_fill_manual(limits = my_levels,values = col) +
      theme_classic() +labs(title = alpha_index[i], x = x.title, y = y.title) +
      theme(
        axis.line = element_line(linewidth = 0.3,lineend = "square"),
        axis.ticks = element_line(linewidth = 0.1),
        plot.title = element_text(size = title.size,hjust = 0.5),
        axis.title.y = element_text(size = axis.title.size,vjust = 1.9,hjust = 0.5,family=family),
        legend.title = element_text(size = legend.title.size,family=family),
        legend.text = element_text(size = legend.text.size,family=family),
        axis.text.x = element_text(size = axis.text.size,vjust = 1, hjust = 1,angle = x.angle,family=family),
        axis.text.y = element_text(size = axis.text.size,vjust = 0.5,hjust = 0.5,family=family)) + 
      scale_y_continuous(labels = scales::scientific)
    
    if (length(my_annotations1) > 0) {
      if (length(levels(plotdata$Group)) <= 5) {
        p <- p +
          geom_signif(annotations = my_annotations1,
                      comparisons = my_comparisons1,
                      step_increase = 0.1,vjust = .2,
                      colour = "gray20",tip_length = 0.015)
      } else {
        x_axis_order <- levels(plotdata$Group)
        y_max <- aggregate(abundance ~ Group, data = plotdata, max)
        y_sep <- diff(range(plotdata$abundance[!is.na(plotdata$abundance)])) * 0.05
        y_start_use <- y_max$abundance + y_sep
        treatments <- as.character(y_max$Group)
        means <- aggregate(abundance ~ Group,data = plotdata, mean)$abundance
        pvalue <- matrix(1,nrow = length(treatments), ncol = length(treatments))
        rownames(pvalue) <- colnames(pvalue) <- treatments
        for (c in my_comparisons1) {
          pvalue[c[1], c[2]] <- 0
          pvalue[c[2], c[1]] <- 0
        }
        grps <- orderPvalue(treatments, means, 0.05, pvalue, console = TRUE)
        add_sig_label <- grps[x_axis_order, "groups",  drop = T]
        
        textdf <- data.frame(x = x_axis_order,y = y_start_use,add = add_sig_label,stringsAsFactors = FALSE)
        
        p <- p +
          geom_text(aes(x = x, y = y, label = add),data = textdf, inherit.aes = FALSE)
      }
    }
    
    ggplotsave(plot=p,
               mapname =paste0(mapname,"-",j),
               savepath = savepath,
               imagetype = imagetype,
               height=height,
               width=width,
               family=family,
               dpi=dpi)
    
  }}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_common_boxplot2 <- map_autodraw$new(moudle = boxplot_new,row.names = 1)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "矩阵文件",
                      required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = "Boxplot", help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf","html"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  parser$add_argument("-cf","--classfile",default = "classtype.xlsx", help = "分组颜色模板")
  
  # 此图参数
  parser$add_argument("-xt","--x.title", default = "", type = "character", help="图片x轴标签")
  parser$add_argument("-yt","--y.title", default = "Intensity", type = "character", help="图片y轴标签")
  parser$add_argument("-go","--grouporder", default = "NA", type = "character", help="分组绘图顺序")
  parser$add_argument("-xa","--x.angle", default = 45, type = "double", help="图片x轴标签角度")
  parser$add_argument("-ts","--title.size", default = 12, type = "double", help="图片标题文字大小")
  parser$add_argument("-ates","--axis.text.size", default = 11, type = "double", help="图片轴标签文字大小")
  parser$add_argument("-atis","--axis.title.size", default = 12, type = "double", help="图片轴标题文字大小")
  parser$add_argument("-ltis","--legend.title.size", default = 11, type = "double", help="图例标签文字大小")
  parser$add_argument("-ltes","--legend.text.size", default = 10, type = "double", help="图例标题文字大小")
  
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  result <- do.call(what = map_common_boxplot2, args = args)
  
}

#' @export
map_common_boxplot2 <- map_autodraw$new(moudle = boxplot_new,row.names = 1)$draw
