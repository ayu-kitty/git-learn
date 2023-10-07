#!/opt/conda/bin/Rscript

#' 棒棒糖图
#'
#' @param data 输入数据
#' @param mapname 图片保存名称
#' @param savepath 图片保存路径
#' @param height 图片保存高度
#' @param width 图片保存宽度
#' @param num 代谢物数量
#' @param point_shape 顶端形状
#' @param x_title_size x轴标题文字大小
#' @param y_title_size y轴标题文字大小
#' @param x_text_size x轴标签文字大小
#' @param y_text_size y轴标签文字大小
#' @param title_size 标题文字大小
#' @param legend_text_size 图例文字大小
#' @param legend_title_size 图例标题大小
#' @param vline_color o值处竖线颜色
#' @param updowncolor 上下调颜色
#' @param title_color 标题颜色
#' @param x_title_color x轴标题颜色
#' @param y_title_color y轴标题颜色
#' @param x_text_color x轴标签颜色
#' @param y_text_color y轴标签颜色
#' @param legend_text_color 图例文字颜色
#' @param legend_title_color 标题颜色
#' @param plot.margin 图片边距
#' @param imagetype 图片保存格式
#' @param title.x x轴标题
#' @param title.y y轴标题
#' @param title.size size图例标题
#' @param title.colour colour图例标题
#' @param plottitle 标题
#' @param ... 见[ggplotsave()]
#' 
#' @export
lolipopmap <- function(data,
                       mapname = "lolipop",
                       savepath = "./lollipopchart/",
                       height = 9,
                       width = 15,
                       num = 10,
                       point_shape = 16,
                       x_title_size = 16,
                       y_title_size = 16,
                       x_text_size = 14,
                       y_text_size = 14,
                       title_size = 16,
                       legend_text_size = 11,
                       legend_title_size = 12,
                       vline_color = "gray",
                       updowncolor = SelectColors(palette = "pinkblue2",n = 2),
                       title_color = "black",
                       x_title_color = "black",
                       y_title_color = "black",
                       x_text_color = "black",
                       y_text_color = "black",
                       legend_text_color = "black",
                       legend_title_color = "black",
                       plot.margin = c(2,2,2,2),
                       imagetype = c("jpg","pdf"),
                       title.x = "Log2(Fold Change)",
                       title.y = "Metabolites",
                       title.size = "VIP",
                       title.colour = "black",
                       aspect.ratio = 16 / 10,
                       family = "sans",
                       plottitle = "",
                       ...){
  
  suppressMessages(library("stringr"))
  suppressMessages(library("ggplot2"))
  if("diff_filter" %in% class(data)){
    data <- get_diffana_data_2(data = data)
  }
  
  drawdata <- data
  
  if("VIP" %in% colnames(drawdata) & "log2FoldChange" %in% colnames(drawdata)){
    
    drawdata$`UPDOWN` <- ifelse(drawdata$`log2FoldChange`>0,"up-regulated","down-regulated")
    drawdata <- drawdata[,c("Metabolites","UPDOWN","VIP","log2FoldChange","p-value")]
    colnames(drawdata) <- c("Metabolites","UPDOWN","VIP","log2(FC)","P-value")
    
    resultdata <- data.frame()
    if("up-regulated" %in% drawdata$UPDOWN){
      updata <- drawdata[drawdata$UPDOWN=="up-regulated",]
      updata <- updata[order(-updata$`VIP`),]
      if(nrow(updata) < num){
        updata <- updata
      }else{
        updata <- head(updata,n=num)
      }
      resultdata <- rbind(resultdata,updata)
      maxlimit <- max(updata$`log2(FC)`)*1.5
    }else{
      maxlimit <- 0
    }
    
    if("down-regulated" %in% drawdata$UPDOWN){
      downdata <- drawdata[drawdata$UPDOWN=="down-regulated",]
      downdata <- downdata[order(downdata$`VIP`),]
      if(nrow(downdata) < num){
        downdata <- downdata
      }else{
        downdata <- tail(downdata,n=num)
      }
      resultdata <- rbind(resultdata,downdata)
      minlimit <- min(downdata$`log2(FC)`)*1.5
    }else{
      minlimit <- 0
    }
    maxlimit <- max(abs(c(maxlimit,minlimit)))
    minlimit <- -max(abs(c(maxlimit,minlimit)))
    
    for(k in 1:nrow(resultdata)){
      resultdata$Metabolites[k] <- unlist(strsplit(split = ";\n",x = resultdata$Metabolites[k]))[1]
      resultdata$Metabolites[k] <- unlist(strsplit(split = "; ",x = resultdata$Metabolites[k]))[1]
      if(str_length(resultdata$Metabolites[k])>35){
        newname <- paste0(substring(resultdata$Metabolites[k],1,35),"...")
        i <- 2
        while (newname %in% resultdata$Metabolites) {
          newname <- paste0(substring(resultdata$Metabolites[k],1,35),"...-",i)
          i <- i+1
        }
        resultdata$Metabolites[k] <- newname 
      }
    }
    
    lolipop <- resultdata
    for(j in 1:nrow(lolipop)){
#      if (resultdata[j,"P-value"]<=0.0001){
#        lolipop[j,"P-value"] <- '****'
#      }else if (resultdata[j,"P-value"]>0.0001&resultdata[j,"P-value"]<=0.001){
#     只显示*;**;***3级别
      if (resultdata[j,"P-value"]<=0.001){
        lolipop[j,"P-value"] <- '***'
      }else if (resultdata[j,"P-value"]>0.001&resultdata[j,"P-value"]<=0.01){
        lolipop[j,"P-value"] <- '**'
      }else if (resultdata[j,"P-value"]>0.01&resultdata[j,"P-value"]<=0.05){
        lolipop[j,"P-value"] <- '*'
      }else{lolipop[j,"P-value"] <- ''}
    }
    
    lolipop$Metabolites <- factor(lolipop$Metabolites,level=lolipop[,1])
    
    pp <- ggplot(lolipop,aes(x=`log2(FC)`,y=Metabolites))+
      geom_segment(aes(x=0,xend=`log2(FC)`,y=Metabolites,yend=Metabolites,colour=UPDOWN),linewidth=2)+
      geom_vline(xintercept=0,linewidth=0.5,colour=vline_color)+
      geom_point(aes(size=VIP,colour=UPDOWN),shape=point_shape)+
      #增加*颜色显示,大小调整
      geom_text(aes(label=`P-value`),hjust = ifelse(lolipop$UPDOWN=="down-regulated", 1.5,-0.5),
                vjust = 0.8,size = 8,color = ifelse(lolipop$UPDOWN=="down-regulated", updowncolor[1],updowncolor[2]))+
      labs(x=title.x,y=title.y,size=title.size,colour="",title=plottitle)+
      scale_color_manual(values = updowncolor )+ 
      scale_y_discrete(expand = expansion(mult=0.05))+
      scale_x_continuous(limit = c(minlimit,maxlimit))+
      theme_bw()+
      guides(size=guide_legend(order = 1))+
      theme(panel.grid = element_blank(),
            legend.title = element_text(size=legend_title_size,color=legend_title_color,angle=0),
            legend.text = element_text(size=legend_text_size,color=legend_text_color,angle=0),
            axis.text.x = element_text(size=x_text_size,color=x_text_color,angle=0),
            axis.text.y = element_text(size=y_text_size,color=y_text_color,angle=0),
            axis.title.y = element_text(size=y_title_size,color=y_title_color,angle=90,vjust=1.5),
            axis.title.x = element_text(size=x_title_size,color=x_title_color,angle=0,vjust=0.5),
            plot.title = element_text(size=title_size,color=title_color,angle=0,hjust=0.5),
            plot.margin = unit(plot.margin, "cm"),
            aspect.ratio = aspect.ratio)
    
    lmbio::ggplotsave(plot = pp,
               savepath = savepath,
               mapname = mapname,
               imagetype = imagetype,
               height = height,
               width = width,
               family = family,
               ...)
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_common_lolipopmap <- map_autodraw$new(moudle = lolipopmap)$draw
  
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
  # parser$add_argument("-ps","--point_shape",default =16, help = "顶端形状")
  # parser$add_argument("-xtis","--x_title_size",default = 14,help = "x轴标题文字大小",dest = "x_title_size")
  # parser$add_argument("-ytis","--y_title_size",default = 14,help = "y轴标题文字大小")
  # parser$add_argument("-xtes","--x_text_size",default = 12,help = "x轴标签文字大小")
  # parser$add_argument("-ytes","--y_text_size",default = 12,help = "y轴标签文字大小")
  # parser$add_argument("-ts","-title_size",default = 16,help = "标题文字大小")
  # parser$add_argument("-ltes","--legend_text_size",default = 11, help = "图例文字大小")
  # parser$add_argument("-ltis","--legend_title_size",default = 12, help = "图例标题大小")
  # parser$add_argument("-vc","--vline_color",default = "gray", help = "o值处竖线颜色")
  # parser$add_argument("-uc","--updowncolor",default = c("#391bfb","#fd1a08"), help = "上下调颜色")
  # parser$add_argument("-tc","--title_color",default = "black", help = "标题颜色")
  # parser$add_argument("-xtic","--x_title_color",default = "black", help = "x轴标题颜色")
  # parser$add_argument("-ytic","--y_title_color",default = "black", help = "y轴标题颜色")
  # parser$add_argument("-xtec","--x_text_color",default = "black", help = "x轴标签颜色")
  # parser$add_argument("-ytec","--y_text_color",default = "black", help = "y轴标签颜色")
  # parser$add_argument("-ltec","--legend_text_color",default ="black", help = "图例文字颜色")
  # parser$add_argument("-ltic","--legend_title_color",default ="black", help = "图例标题颜色")
  # parser$add_argument("-pm","--plot.margin",default =c(2,2,2,2), help = "绘图边距")
  # parser$add_argument("-tx","--title.x",default ="Log2(Fold Change)", help = "x轴标题")
  # parser$add_argument("-ty","--title.y",default ="Metabolites", help = "y轴标题")
  # parser$add_argument("-tis","--title.size",default ="VIP", help = "size图例标题")
  # parser$add_argument("-tic","--title.colour",default ="", help = "colour图例标题")
   parser$add_argument("-n","--num",default = 10, help = "上调或者下调代谢物数量")
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  result <- do.call(map_common_lolipopmap,args = args)
  
  warnings()
}

map_common_lolipopmap <- map_autodraw$new(moudle = lolipopmap)$draw
