#!/opt/conda/bin/Rscript

#' 分类气泡图
#'
#' @param data 数据
#' @param mapname 图片保存名称
#' @param savepath 图片保存路径
#' @param height 图片保存高度
#' @param width 图片保存高度
#' @param x_title_size x轴标题文字大小
#' @param y_title_size y轴标题文字大小
#' @param y_text_size y轴标签文字大小
#' @param title_size 标题文字大小
#' @param legend_text_size 图例文字大小
#' @param legend_title_size 图例标题大小
#' @param family 字体,默认为Aria
#' @param title_color  绘图标题颜色
#' @param x_title_color 轴标题颜色
#' @param y_title_color y轴标题颜色
#' @param y_text_color y轴标签颜色
#' @param legend_text_color 图例文字颜色
#' @param legend_title_color 图例标题颜色
#' @param plot.margin 绘图边距
#' @param imagetype 图片保存格式
#' @param title.x x轴标题
#' @param title.y y轴标题
#' @param title.size size图例标题
#' @param color_class Class颜色
#' @param sub_colorlist SubClass颜色
#' 
#' @export
classfication_bubbleplot <- function(data,
                                     mapname="bubbleplot",
                                     savepath="./脂质分类气泡图/",
                                     height=7,
                                     width=10,
                                     family="sans",
                                     x_title_size=14,
                                     y_title_size=14,
                                     y_text_size=13,
                                     title_size=16,
                                     legend_text_size=11,
                                     legend_title_size=12,
                                     title_color="black",
                                     x_title_color="black",
                                     y_title_color="black",
                                     y_text_color="black",
                                     legend_text_color="black",
                                     legend_title_color="black",
                                     plot.margin=c(2,2,2,2),
                                     imagetype=c("png","pdf","jpg"),
                                     title.x="Lipid Classfication",
                                     title.y="Log2(Fold Change)",
                                     title.size="P-value",
                                     color_class=c("#880808","#DE3163","#0F52BA","#097969","#FFBF00","#F28C28","#A52A2A"),
                                     sub_colorlist=list(c("#A52A2A","#C41E3A","#800020","#DC143C","#D2042D",
                                                          "#E97451","#AA4A44","#D22B2B","#D70040","#EE4B2B",
                                                          "#811331","#F88379","#6E260E","#DE3163","#8B0000",
                                                          "#7B1818","#9A2A2A","#C04000","#800000","#722F37",
                                                          "#770737","#FF3131","#4A0404","#FAA0A0","#EC5800",
                                                          "#E35335","#E30B5C","#FF0000","#A52A2A","#A42A04",
                                                          "#C21E56","#E0115F","#FA8072","#FF2400","#E34234"),
                                                        c("#9F2B68","#C21E56","#E0115F","#E30B5C","#F33A6A",
                                                          "#FF69B4","#E37383","#FFB6C1","#AA336A","#FF00FF",
                                                          "#FF10F0","#DA70D6","#702963","#FAA0A0","#FFC0CB",
                                                          "#F89880","#673147","#A95C68","#800080","#FF00FF",
                                                          "#953553","#E0BFB8","#811331","#F8C8DC","#F88379",
                                                          "#FA8072","#E0B0FF","#D8BFD8","#BF40BF","#E34234",
                                                          "#CBC3E3","#CF9FFF","#AA98A9","#FFF5EE","#F3CFC6"),
                                                        c("#00008B","#0000FF","#0047AB","#6495ED","#89CFF0",
                                                          "#1F51FF","#0096FF","#4169E1","#1434A4","#0F52BA",
                                                          "#6F8FAF","#7DF9FF","#6082B6","#00FFFF","#000080",
                                                          "#3F00FF","#5D3FD3","#ADD8E6","#191970","#0437F2",
                                                          "#A7C7E7","#CCCCFF","#B6D0E2","#96DED1","#088F8F",
                                                          "#9FE2BF","#87CEEB","#4682B4","#7393B3","#5F9EA0",
                                                          "#008080","#40E0D0","#40B5AD","#0818A8","#F0FFFF"),
                                                        c("#355E3B","#4F7942","#008000","#2E8B57","#228B22",
                                                          "#AFE1AF","#50C878","#32CD32","#4CBB17","#0BDA51",
                                                          "#7CFC00","#00A36C","#2AAA8A","#088F8F","#00FFFF",
                                                          "#90EE90","#478778","#454B1B","#AAFF00","#5F8575",
                                                          "#98FB98","#8A9A5B","#0FFF50","#40E0D0","#808000",
                                                          "#C1E1C1","#B4C424","#93C572","#8A9A5B","#7FFFD4",
                                                          "#9FE2BF","#009E60","#00FF7F","#008080","#ECFFDC"),
                                                        c("#FDDA0D","#FFEA00","#FCF55F","#FFDB58","#E4D00A",
                                                          "#DFFF00","#F4C430","#FFFF00","#FFAA33","#FFFF8F",
                                                          "#C2B280","#EEDC82","#E49B0F","#FFD700","#FFC000",
                                                          "#DAA520","#F8DE7E","#F0E68C","#FAFA33","#FFF8DC",
                                                          "#FBEC5D","#F4BB44","#FADA5E","#FFDEAD","#FAD5A5",
                                                          "#FFFAA0","#FFE5B4","#C9CC3F","#F3E5AB","#E1C16E",
                                                          "#C4B454","#F5DEB3","#B4C424","#8B8000","#FFFDD0"),
                                                        c("#FF7518","#FF5F15","#CC5500","#EC5800","#FFAC1C",
                                                          "#FFBF00","#F08000","#F4BB44","#E3963E","#E49B0F",
                                                          "#DAA520","#FFD580","#FF5F1F","#B87333","#FFAA33",
                                                          "#CC7722","#FFA500","#FAC898","#F89880","#FAD5A5",
                                                          "#FF7518","#FF4433","#CD7F32","#FA8072","#FA5F55",
                                                          "#FFC000","#FAFA33","#FCF55F","#FFEA00","#E4D00A",
                                                          "#E3735E","#FFD700","#F4C430","#8B4000","#D27D2D"),
                                                        c("#834333","#7B3F00","#6E260E","#CD7F32","#954535",
                                                          "#B87333","#814141","#6F4E37","#C19A6B","#80461B",
                                                          "#D27D2D","#800020","#5C4033","#8B4513","#A0522D",
                                                          "#8B0000","#988558","#C2B280","#C19A6B","#913831",
                                                          "#E5AA70","#E5AA70","#966919","#C4A484","#C04000",
                                                          "#800000","#967969","#F2D2BD","#CC7722","#A95C68",
                                                          "#D2B48C","#EADDCA","#E1C16E","#E97451","#DAA06D")),
                                     ...){

  pacman::p_load(gridExtra,ggplot2,dplyr)
  
  getLegend <- function(p){
    g <- ggplotGrob(p) 
    k <- which(g$layout$name=="guide-box") 
    g$grobs[[k]] 
  }
  
  if("diff_filter" %in% class(data)){
    data <- get_diffana_data_2(data = data)
  }
  
  data <- data[,c("Metabolites","p-value","log2FoldChange","Sub Class","Class")]
  colnames(data) <- c("Lipid","P-value","Log2(FC)","Sub Class","Class")
  
  data <- data[order(data$Class,data$`Sub Class`,data$`P-value`),]
  dataclass <- distinct(data[,c("Sub Class","Class")])
  Class <- unique(dataclass$`Class`)
  datatotal <- data.frame()
  for(j in 1:length(Class)){
    dataclass_j=dataclass[dataclass$`Class`==Class[j],]
    for(k in 1:nrow(dataclass_j)){
      dataclass_j$`color`[k]<-sub_colorlist[[j]][k]
    }
    datatotal<-rbind(datatotal,dataclass_j)
  }
  collist <- setNames(datatotal$color, datatotal$`Sub Class`)
  pp <- list()
  
  for(i in 1:length(Class)){
    p_i <- ggplot(data,aes(x=factor(Lipid,level=Lipid),y=`Log2(FC)`,fill=`Sub Class`))+
      geom_point(shape=21,color="black") +
      scale_fill_manual(name=Class[i],breaks=dataclass[dataclass$Class==Class[i],1],values=collist)+
      theme_bw()+geom_hline(yintercept=0,linetype=2,color="gray")+
      guides(fill=guide_legend(override.aes = list(size=3)))+
      theme(panel.grid=element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.title =element_text(size=legend_title_size,angle=0,family=family,color=color_class[i]),
            legend.text =element_text(size=legend_text_size,angle=0,family=family),
            axis.text.y=element_text(size=y_text_size,color=y_text_color,angle=0,family=family),
            axis.title.y=element_text(size=y_title_size,color=y_title_color,angle=90,vjust=1.5,family=family),
            axis.title.x=element_text(size=x_title_size,color=x_title_color,angle=0,vjust=0.5,family=family),
            plot.title=element_text(size=title_size,color=title_color,angle=0,hjust=0.5,family=family),
            plot.margin = unit(plot.margin, "cm"),
            panel.border=element_blank(),
            axis.line = element_line())
    
    pp[[i]] <- getLegend(p_i)
  }
  
  
  p_p <- ggplot(data,aes(x=factor(Lipid,level=Lipid),y=`Log2(FC)`,size=`P-value`))+
    geom_point(color="black") +
    theme_bw()+geom_hline(yintercept=0,linetype=2,color="gray")+labs(size="P-value")+
    theme(panel.grid=element_blank(),
          legend.key.size  = unit(5,"pt"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title =element_text(size=legend_title_size,angle=0,family=family),
          legend.text =element_text(size=legend_text_size,angle=0,family=family),
          axis.text.y=element_text(size=y_text_size,color=y_text_color,angle=0,family=family),
          axis.title.y=element_text(size=y_title_size,color=y_title_color,angle=90,vjust=1.5,family=family),
          axis.title.x=element_text(size=x_title_size,color=x_title_color,angle=0,vjust=0.5,family=family),
          plot.title=element_text(size=title_size,color=title_color,angle=0,hjust=0.5,family=family),
          plot.margin = unit(plot.margin, "cm"),
          panel.border=element_blank(),
          axis.line = element_line())
  
  #拼接图例
  pp[[length(Class)+1]] <- getLegend(p_p)
  #do.call(grid.arrange, c(pp, ncol=1))
  maxWidth1 <- pp[[1]]$widths[2:5]
  for(i in 1:length(pp)) {
    maxWidth <- grid::unit.pmax(maxWidth1, pp[[i]]$widths[2:5])
  }
  for(i in 1:length(pp)) {
    pp[[i]]$widths[2:5] <- maxWidth
  }
  #do.call(grid.arrange, c(pp, ncol=1))
  
  p <- ggplot(data,aes(x=factor(Lipid,level=Lipid),y=`Log2(FC)`,size=`P-value`,fill=`Sub Class`))+
    geom_point(shape=21,color="black") +
    scale_fill_manual(name="",values=collist)+
    theme_bw()+geom_hline(yintercept=0,linetype=2,color="gray")+
    scale_x_discrete(expand=expansion(mult=c(0.05,0.1)))+
    scale_y_continuous(expand=expansion(mult=c(0.05,0.1)))+
    labs(size=title.size,x=title.x,y=title.y,title="")+
    theme(panel.grid=element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title =element_text(size=legend_title_size,angle=0,family=family),
          legend.text =element_text(size=legend_text_size,angle=0,family=family),
          axis.text.y=element_text(size=y_text_size,color=y_text_color,angle=0,family=family),
          axis.title.y=element_text(size=y_title_size,color=y_title_color,angle=90,vjust=1.5,family=family),
          axis.title.x=element_text(size=x_title_size,color=x_title_color,angle=0,vjust=0.5,family=family),
          plot.title=element_text(size=title_size,color=title_color,angle=0,hjust=0.5,family=family),
          plot.margin = unit(plot.margin, "cm"),
          panel.border=element_blank(),
          axis.line = element_line())
  
  #添加图例
  g <- ggplotGrob(p)
  k <- which(g$layout$name=="guide-box")
  g$grobs[[k]] <- do.call(grid.arrange, c(pp, ncol=1))
  
  graphics.off()
  ggplotsave(plot=g,
             savepath = savepath,
             mapname = mapname,
             imagetype = imagetype,
             height = height,
             width = width,
             family = family,
             ...)
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  
  map_lipid_classfication_bubbleplot <- map_autodraw$new(moudle = classfication_bubbleplot)$draw
  
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
  # parser$add_argument("-xtis","--x_title_size",default = 14,help = "x轴标题文字大小")
  # parser$add_argument("-ytis","--y_title_size",default = 14,help = "y轴标题文字大小")
  # parser$add_argument("-ytes","--y_text_size",default = 12,help = "y轴标签文字大小")
  # parser$add_argument("-ts","-title_size",default = 16,help = "标题文字大小")
  # parser$add_argument("-ltes","--legend_text_size",default = 11, help = "图例文字大小")
  # parser$add_argument("-ltis","--legend_title_size",default = 12, help = "图例标题大小")
  # parser$add_argument("-tic","--title_color",default = "black", help = "绘图标题颜色")
  # parser$add_argument("-xtic","--x_title_color",default = "black", help = "x轴标题颜色")
  # parser$add_argument("-ytic","--y_title_color",default = "black", help = "y轴标题颜色")
  # parser$add_argument("-ytec","--y_text_color",default = "black", help = "y轴标签颜色")
  # parser$add_argument("-ltec","--legend_text_color",default ="black",help = "图例文字颜色")
  # parser$add_argument("-ltic","--legend_title_color",default ="black", help = "图例标题颜色")
  # parser$add_argument("-pm","--plot.margin",default =16, help = "绘图边距")
  # parser$add_argument("-tx","--title.x",default ="Lipid species", help = "x轴标题")
  # parser$add_argument("-ty","--title.y",default ="Log2(Fold Change)", help = "y轴标题")
  # parser$add_argument("-tis","--title.size",default ="P-value", help = "size图例标题")
  # parser$add_argument("-cc","--color_class",default =c("#a76283","#c67915","#354e6b","#446a37","#e60012"), help = "Class颜色")
  # parser$add_argument("-sc","--sub_colorlist",default =list(c("#dd7694","#ecb0c1","#f9d3e3"),
  #                                                           c("#db9b34","#fac03d","#fedc5e"),
  #                                                           c("#284852","#007175","#3271ae"),
  #                                                           c("#4c8045","#68945c","#a8bf8f"),
  #                                                           c("#e94829","#ed6d3d","#f29a76")), help = "SubClass颜色")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  result <- do.call(map_lipid_classfication_bubbleplot,args = args)
}

#' @export
map_lipid_classfication_bubbleplot <- map_autodraw$new(moudle = classfication_bubbleplot)$draw
