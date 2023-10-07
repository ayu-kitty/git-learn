#! /opt/conda/bin/Rscript

#' 火山图
#'
#' @param data 数据
#' @param mapname 保存文件名
#' @param imagetype 保存图片格式
#' @param height 图片长度
#' @param width 图片宽度
#' @param x x轴标签
#' @param y y轴标签
#' @param showname 逻辑，是否显示名称
#' @param plot.title 图标题主题
#' @param color 颜色
#' @param namelist 名称列名
#' @param labelname 希望标识的名称
#' @param correspondname 对应列表
#' @param max.overlaps 名称重叠度
#' @param drawFC FC标准
#' @param drawFC_linecolor FC筛选线颜色
#' @param drawFC_linetype FC筛选线类型
#' @param cex 点大小
#' @param aspect.ratio 长宽比
#' @param size vip情况下点大小范围
#' @param plot.margin 图像空白边框
#' @param legend.position 图例位置
#' @param legend.background 图例背景
#' @param VIP VIP筛选标准
#' @param p p筛选标准
#' @param linecolor p筛选线颜色
#' @param linetype p筛选线类型
#' @param adjust 逻辑，是否校正p值
#' @param adjustoutx 超出x轴范围数值调整
#' @param adjustouty 超出y轴范围数值调整
#' @param ... 
#'
#' @export
auto_volplot <- function(data,
                         mapname = "volcano",
                         imagetype = c("png","pdf","html"),
                         width = 7,
                         height = 6,
                         x = min(c(max(abs(c(floor(quantile(data[, "log2(FC)"], probs = c(0.999)) * 1.2),
                                             floor(quantile(data[, "log2(FC)"], probs = c(0.001)) * 1.2), 1))), 12)),
                         y = min(c(max(abs(c(floor(quantile(data[, "-log10p"], probs = c(0.999)) * 2), 4))), 50)),
                         drawFC = NULL, drawFC_linecolor = "black", drawFC_linetype = "dashed",
                         cex = 1,
                         aspect.ratio = 1,
                         size = c(0.5, 3),
                         # color = c("blue","cyan", "grey","pink", "red"),
                         color = SelectColors(palette = "volcanocol",n = 5),
                         plot.margin = ggplot2::unit(c(0.3, 0.3, 0.3, 0.3), "in"),
                         plot.title = ggplot2::element_text(hjust = 0.5),
                         legend.position = "right",
                         legend.background = ggplot2::element_rect(),
                         VIP = NULL,
                         p = 0.05,
                         linecolor = "black",
                         linetype = "dashed",
                         adjust = F,
                         adjustoutx = T,
                         adjustouty = T,
                         showname = F,
                         namenum = 5,
                         namelist = NULL,
                         labelname = NULL,
                         correspondname = NULL,
                         max.overlaps = 200,
                         ...) {
  
  suppressMessages(library("ggplot2"))
  suppressMessages(library("ggrepel"))
  
  options(warn = -1)
  
  vol <- data
  
  # 调整出界点
  if (adjustouty) {vol[vol[, "-log10p"] > y, "-log10p"] <- y}
  if (adjustoutx) {vol[vol[, "log2(FC)"] > x, "log2(FC)"] <- x
  vol[vol[, "log2(FC)"] < -x, "log2(FC)"] <- -x}
  
  #### 火山图作图代码####
  vol[, "class"] <- "not-significant"
  
  if (is.null(drawFC)) {
    if(is.null(p)){
      vol$class[vol[, "log2(FC)"] > 0] <- "up-regulated"
      vol$class[vol[, "log2(FC)"] < 0] <- "down-regulated"
      upfiltertxt <- paste0("log2(FC) > 0","\n","Up")
      downfiltertxt <- paste0("log2(FC) < 0","\n","Down")
    }else{
      vol$class[vol[, "-log10p"] > -log10(p[1]) & vol[, "log2(FC)"] > 0] <- "up-regulated"
      vol$class[vol[, "-log10p"] > -log10(p[1]) & vol[, "log2(FC)"] < 0] <- "down-regulated"
      if(adjust){
        upfiltertxt <- paste0("q < ",p[1],"\n","log2(FC) > 0","\n","Up")
        downfiltertxt <- paste0("q < ",p[1],"\n","log2(FC) < 0","\n","Down")
      }else{
        upfiltertxt <- paste0("p < ",p[1],"\n","log2(FC) > 0","\n","Up")
        downfiltertxt <- paste0("p < ",p[1],"\n","log2(FC) < 0","\n","Down")
      }
    }
    
    breaks <- c("down-regulated","not-significant","up-regulated")
    labels <- c("Significant Down","Non-significant","Significant Up")
  } else {
    if(is.null(p)){
      vol$class[vol[, "log2(FC)"] > 0] <- "up-regulated-notsign"
      vol$class[vol[, "log2(FC)"] < 0] <- "down-regulated-notsign"
      vol$class[vol[, "log2(FC)"] > log2(drawFC)] <- "up-regulated"
      vol$class[vol[, "log2(FC)"] < (-log2(drawFC))] <- "down-regulated"
      upfiltertxt <- paste0("log2(FC) > ",round(log2(drawFC),3),"\n","Up")
      downfiltertxt <- paste0("log2(FC) < -",round(log2(drawFC),3),"\n","Down")
    }else{
      vol$class[vol[, "-log10p"] > -log10(p[1]) & vol[, "log2(FC)"] > 0] <- "up-regulated-notsign"
      vol$class[vol[, "-log10p"] > -log10(p[1]) & vol[, "log2(FC)"] < 0] <- "down-regulated-notsign"
      vol$class[vol[, "-log10p"] > -log10(p[1]) & vol[, "log2(FC)"] > log2(drawFC)] <- "up-regulated"
      vol$class[vol[, "-log10p"] > -log10(p[1]) & vol[, "log2(FC)"] < (-log2(drawFC))] <- "down-regulated"
      if(adjust){
        upfiltertxt <- paste0("q < ",p[1],"\n","log2(FC) > ",round(log2(drawFC),3),"\n","Up")
        downfiltertxt <- paste0("q < ",p[1],"\n","log2(FC) < -",round(log2(drawFC),3),"\n","Down")
      }else{
        upfiltertxt <- paste0("p < ",p[1],"\n","log2(FC) > ",round(log2(drawFC),3),"\n","Up")
        downfiltertxt <- paste0("p < ",p[1],"\n","log2(FC) < -",round(log2(drawFC),3),"\n","Down")
      }
    }
    
    if(color[2] == color[3] & color[4] == color[3]){
      breaks <- c("down-regulated","not-significant","up-regulated")
      labels <- c("Significant Down","Non-significant","Significant Up")
    }else{
      breaks <- c("down-regulated","down-regulated-notsign","not-significant","up-regulated-notsign","up-regulated")
      labels <- c("Significant Down", "Down","Non-significant","Up","Significant Up")
    }
    
  }
  
  if(length(VIP)==1){
    if(VIP == 0){
      VIP <- NULL
    }
  }
  if (is.null(VIP)) {
    
    upsum <- sum(vol$class == "up-regulated")
    downsum <- sum(vol$class == "down-regulated")
    
    pp <- ggplot(data = vol,
                 mapping = aes(x = `log2(FC)`, y = `-log10p`, color = class,
                               text = paste0("ID:", row.names(vol)))) +
      geom_point(size = cex) +
      scale_x_continuous(limits = c(-x, x)) +
      scale_y_continuous(limits = c(0, y)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_rect(fill=NA,color="black", linewidth=1),
            plot.title = plot.title,
            aspect.ratio = aspect.ratio,
            plot.margin = plot.margin,
            legend.position = legend.position,
            legend.background = legend.background,
            line = element_line(colour = "black"),
            legend.margin = margin(t = 0, r = 6, b = 6, l = 2, unit = "pt")) +
      #火山图新增文字部分
      annotate("text",label=paste0("Sig : ",upsum,"\n",
                                   upfiltertxt),
               x=x-0.2,y=y-y/10,size=3,color="black",hjust=1)+
      annotate("text",label=paste0("Sig : ",downsum,"\n",
                                   downfiltertxt),
               x=-x+0.2,y=y-y/10,size=3,color="black",hjust=0)+
      guides(color = guide_legend(title = NULL)) +
      scale_color_manual(values = c("down-regulated" = color[1],
                                    "down-regulated-notsign" = color[2],
                                    "not-significant" = color[3],
                                    "up-regulated-notsign" = color[4],
                                    "up-regulated" = color[5]),
                         breaks = breaks,
                         labels = labels)
    
  } else {
    maxn <- max(vol[, "VIP"],na.rm = T)
    
    upfiltertxt <- paste0("VIP > ",VIP,"\n",upfiltertxt)
    downfiltertxt <- paste0("VIP > ",VIP,"\n",downfiltertxt)
    
    if (maxn < 1) { lab <- c(0.3, 0.6, 0.9)
    } else if (maxn < 3) { lab <- c(0.5, 1, 2)
    } else if (maxn < 5) { lab <- c(0.5, 1, 2, 4)
    } else if (maxn < 8) { lab <- c(1, 3, 5, 7)
    } else if (maxn < 10) { lab <- c(1, 3, 5, 7, 9)
    } else if (maxn < 20) { lab <- c(1, 6, 11, 16)
    } else if (maxn < 30) { lab <- c(1, 9, 17, 25)
    } else if (maxn < 50) { lab <- c(1, 12, 24, 35)
    } else { lab <- c(1, 15, 30, 45) }
    
    vol$class[vol$VIP < VIP] <- "not-significant"
    vol <- vol[order(vol$VIP), ]
    
    upsum <- sum(vol$class == "up-regulated")
    downsum <- sum(vol$class == "down-regulated")
    
    pp <- ggplot(data = vol,
                 mapping = aes(x = `log2(FC)`, y = `-log10p`,
                               color = class, size = VIP,
                               text = paste0("ID:", row.names(vol)))) +
      geom_point() +
      scale_size_continuous(breaks = lab, range = size) +
      scale_x_continuous(limits = c(-x, x)) +
      scale_y_continuous(limits = c(0, y)) +
      theme_bw() +
      #火山图新增文字部分
      annotate("text",label=paste0("Sig : ",upsum,"\n",
                                   upfiltertxt),
               x=x-0.2,y=y-y/10,size=3,color="black",hjust=1)+
      annotate("text",label=paste0("Sig : ",downsum,"\n",
                                   downfiltertxt),
               x=-x+0.2,y=y-y/10,size=3,color="black",hjust=0)+
      theme(panel.grid = element_blank(),
            panel.border = element_rect(fill=NA,color="black", linewidth=1),
            plot.title = plot.title,
            aspect.ratio = aspect.ratio,
            plot.margin = plot.margin,
            legend.position = legend.position,
            legend.background = legend.background,
            line = element_line(colour = "black"),
            legend.margin = margin(t = 0, r = 6, b = 6, l = 2, unit = "pt")) +
      guides(color = guide_legend(title = NULL)) +
      scale_color_manual(values = c("down-regulated" = color[1],
                                    "down-regulated-notsign" = color[2],
                                    "not-significant" = color[3],
                                    "up-regulated-notsign" = color[4],
                                    "up-regulated" = color[5]),
                         breaks = breaks,
                         labels = labels)+
      guides(size = guide_legend(order = 1))
  }
  
  if (adjust) {
    pp <- pp + labs(x = "log2(FC)",
                    y = "-log10(q-value)",
                    title = "Volcano Plot")
  } else {
    pp <- pp + labs(x = "log2(FC)",
                    y = "-log10(P-value)",
                    title = "Volcano Plot")
  }
  
  if (is.null(p)) {
  } else {
    for (i in 1:length(p)) {
      pp <- pp + geom_hline(yintercept = c(-log10(p[i])),
                            linetype = linetype,
                            size = 0.5,
                            colour = linecolor) +
        annotate("text",
                 label = paste0(ifelse(adjust,"q-value","p-value"),"=", p[i]),
                 x = ifelse(adjust,0.8,0.8)*x, y = c(-log10(p[i])) - y * 0.02,
                 size = 4, colour = linecolor)
    }
  }
  
  
  if (is.null(drawFC)) {
  } else {
    pp <- pp + geom_vline(xintercept = c(log2(drawFC)),
                          linetype = drawFC_linetype,
                          size = 0.5,
                          colour = drawFC_linecolor) +
      geom_vline(xintercept = c(-log2(drawFC)),
                 linetype = drawFC_linetype,
                 size = 0.5,
                 colour = drawFC_linecolor)# +
    # annotate("text",
    #          label = paste0("FC=", drawFC),
    #          x = log2(drawFC) + 0.05 * x,
    #          y = y * 0.95,
    #          size = 4,
    #          colour = drawFC_linecolor,
    #          hjust = 0) +
    # annotate("text",
    #          label = paste0("FC=", "1/", drawFC),
    #          x = -log2(drawFC) - 0.05 * x,
    #          y = y * 0.95,
    #          size = 4,
    #          colour = drawFC_linecolor,
    #          hjust = 1)
  }
  
  if (showname) {
    vol[, "name"] <- NA
    
    if (is.null(namelist)) {
      vol[, "name2"] <- vol[, 1]
    } else {
      vol[, "name2"] <- vol[, namelist] 
    }
    
    if (is.null(labelname)) {
    } else {
      vol[vol[, "name2"] %in% labelname, "name"] <- vol[vol[, "name2"] %in% labelname, "name2"]
      if (!is.null(correspondname)) {vol[, "name"] <- whto(correspondname, vol[, "name"])}
    }
    
    if(!is.null(namenum)){
      if(namenum>0){
        vol2 <- vol[vol$class == "down-regulated",]
        if(dim(vol2)[1] > 0){
          if(dim(vol2)[1] < namenum){
            namenum2 <- dim(vol2)[1]
          }else{
            namenum2 <- namenum
          }
          vol2 <- vol2[order(vol2$`log2(FC)`),]
          vol[row.names(vol2)[1:namenum2],"name"] <- vol[row.names(vol2)[1:namenum2],"name2"]
        }
        vol2 <- vol[vol$class == "up-regulated",]
        if(dim(vol2)[1] > 0){
          if(dim(vol2)[1] < namenum){
            namenum2 <- dim(vol2)[1]
          }else{
            namenum2 <- namenum
          }
          vol2 <- vol2[order(vol2$`log2(FC)`,decreasing = T),]
          vol[row.names(vol2)[1:namenum2],"name"] <- vol[row.names(vol2)[1:namenum2],"name2"]
        }
      }
    }
    
    pp <- pp + ggnewscale::new_scale_color()+
      geom_point(data = vol[!is.na(vol$name),], 
                 mapping =  aes(x = `log2(FC)`, y = `-log10p`,label = name,color = class),
                 inherit.aes = F,show.legend = F,size = cex)+
      ggrepel::geom_label_repel(data = vol[!is.na(vol$name),], 
                                mapping = aes(x = `log2(FC)`, y = `-log10p`,label = name,color = class),
                                size = 2.5,
                                box.padding = 0.5,
                                min.segment.length = 0.1,
                                max.overlaps = max.overlaps,
                                inherit.aes = F,
                                na.rm = T,
                                show.legend = F)+
      scale_color_manual(values = c("down-regulated" = addblack_color(color[1]),
                                    "down-regulated-notsign" = addblack_color(color[2]),
                                    "not-significant" = addblack_color(color[3]),
                                    "up-regulated-notsign" = addblack_color(color[4]),
                                    "up-regulated" = addblack_color(color[5])),
                         breaks = breaks,
                         labels = labels)
  }
  
  args <- ggplotsave(plot = pp,
                     mapname = mapname,
                     width = width,
                     height = height,
                     imagetype = imagetype,
                     ...)
  
  return(args)
}


#' 火山图
#'
#' @param data 数据
#' @param mapname 保存文件名
#' @param ... 见[auto_volplot()]
#'
#' @export
auto_volcano <- function(data,
                         mapname ="Volcano",
                         adjust = F,
                         ...) {
  vol <- data
  
  if(adjust){
    pvaluename <- c("q-value","adj.P-value","FDR")
  }else{
    pvaluename <- c("P-value","p-value")
  }
  
  for ( i in 1:length(pvaluename)) {
    if(pvaluename[i] %in% colnames(vol)){
      vol[, "-log10p"] <- (-log10(vol[, pvaluename[i]]))
      break
    }
    
    if(i == length(pvaluename)){
      stop("未找到pvalue相关列")
    }
  }
  
  log2fcname <- c("log2(FC)","log2FoldChange")
  
  for ( i in 1:length(log2fcname)) {
    
    if(log2fcname[i] %in% colnames(vol)){
      vol[, "log2(FC)"] <- vol[, log2fcname[i]]
      break
    }
    
    if(i == length(log2fcname)){
      fcname <- c("FC","FoldChange")
      for ( j in 1:length(fcname)) {
        
        if(log2fcname[i] %in% colnames(vol)){
          
          vol[, "log2(FC)"] <- log2(vol[, fcname[i]])
          
        }
        if(i == length(log2fcname)){
          stop("未找到fc相关列")
        }
      }
    }
  }
  
  vol <- vol[!is.na(vol[, "-log10p"]), ]
  vol <- vol[!is.na(vol[, "log2(FC)"]), ]
  vol <- vol[!vol[, "log2(FC)"] == 0, ]
  
  returndata <- auto_volplot(data = vol,
                             mapname = mapname,
                             adjust = adjust,
                             ...)
  
  # print("auto_volcano运行完成")
  
  return(returndata)
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_common_volcano <- map_autodraw$new(auto_volcano)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "火山图矩阵文件",required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("png","pdf","html"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf","html"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-pf","--pfilter", default =  0.05,type = "double",
                      help = "pvalue的筛选标准",dest = "p")
  parser$add_argument("-ff","--fcfilter", default =  0,type = "double",
                      help = "fc的筛选标准",dest = "drawFC")
  parser$add_argument("-vf","--vipfilter", default =  0,type = "double",
                      help = "vip的筛选标准",dest = "VIP")
  parser$add_argument("-c","--color", default = "volcanocol", 
                      help = "散点颜色，默认蓝、浅蓝、灰、浅红、红")
  parser$add_argument("-sn","--showname", default =  "F", help = "是否显示名称")
  parser$add_argument("-nn","--namenum", default =  0,type = "integer",help = "显示名称数量")
  parser$add_argument("-nl","--namelist", default =  "", help = "显示名称使用的列")
  parser$add_argument("-ln","--labelname",default = "", help = "需标记名字",nargs = "+")
  parser$add_argument("-ad","--adjust",default = F, help = "是否使用校正的p值",action='store_true')
  
  args <- parser$parse_args()
  
  args$color <- SelectColors(palette = args$color,n = 5)
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  if(args$drawFC == 0){args$drawFC <- NULL}
  if(args$VIP == 0){args$VIP <- NULL}
  if(args$showname == "T"){args$showname <- T}else{args$showname <- F}
  if(args$namelist == ""){args$namelist <- NULL}
  
  if(grepl(pattern = ".xlsx$",args$labelname[1]) |
     grepl(pattern = ".xls$",args$labelname[1]) |
     grepl(pattern = ".txt$",args$labelname[1])|
     grepl(pattern = ".csv$",args$labelname[1])){
    if(file.exists(args$labelname[1])){
      labelname <- readdata(args$labelname[1])
      labelname <- c(colnames(labelname)[1],labelname[,1])
      args$labelname <- labelname
      if(labelname[1] == "Gene Name"){args$namelist <- "Gene Name"}
      if(labelname[1] == "Accession"){args$namelist <- "Accession"}
      if(labelname[1] == "Metabolites"){args$namelist <- "Metabolites"}
      # base::print(args$labelname)
    }
  }
  
  result <- do.call(what = map_common_volcano,args = args) 
  
}

#' 根据文件进行火山图可视化
#' 
#' @export
map_common_volcano <- map_autodraw$new(auto_volcano)$draw
