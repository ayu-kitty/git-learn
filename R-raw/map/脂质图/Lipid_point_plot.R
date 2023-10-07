#!/opt/conda/bin/Rscript

#' @export
Lipid_point_plot <- function(data,
                             plot.caption_size = 19,
                             axis.text_size = 14,
                             axis.title_size = 16,
                             legend.title_size = 15,
                             legend.text_size = 10,
                             mapname = "绘图数据-point.xlsx",
                             imagetype = c("jpg", "pdf"),
                             savepath = "./",
                             family = "sans",
                             width = 12,
                             height = 8,
                             dpi = 300,
                             wbsave = T,
                             savename = paste0(mapname,".xlsx"),
                             ...){
  suppressMessages(library(reshape2))
  suppressMessages(library(ggplot2))
  
  if("diff_filter" %in% class(data)){
    mat <- data
    plotdata <- merge(mat$args$info[,c("Sub Class","Total Carbon","Total Unsaturation"),drop=F],
                      mat$args$data,by = 0,all = F)[,-1]
    data_C <- aggregate(.~`Sub Class`+`Total Carbon`+`Total Unsaturation`,plotdata,sum) 
    
    diffdata <- fc_ana(data = data_C[,-1:-2,drop = F],
                       class = mat$args$singleclass,
                       compare = mat$args$compare,
                       log = mat$args$log)
    
    diffdata2 <- p_ana(data = data_C[,-1:-2,drop = F],
                       class = mat$args$singleclass,
                       paired = mat$args$paired,
                       p_method = mat$args$p_method,
                       p_adjust_methods = mat$args$p_adjust_methods)
    diffdata <- merge(diffdata,diffdata2,by = 0,sort = F,all = T)
    row.names(diffdata) <- diffdata[,1]
    diffdata <- diffdata[,-1]
    data_C <- merge(data_C,diffdata,by = 0,sort = F,all = T)[,-1]
    result <- data_C
  }else{
    result <- data
  }
  
  Class <- unique(result$`Sub Class`)
  Unsaturation <- unique(result$`Total Unsaturation`)
  
  for(k in 1:length(Class)){
    
    class_data <- result[result$`Sub Class`==Class[k],]
    class_data$`Total Carbon` <- as.character(class_data$`Total Carbon`)
    class_data$`Total Carbon` <- factor(class_data$`Total Carbon`,levels = unique(class_data$`Total Carbon`))
    class_data$`Total Unsaturation` <- as.character(class_data$`Total Unsaturation`)
    class_data$`Total Unsaturation` <- factor(class_data$`Total Unsaturation`,levels = unique(class_data$`Total Unsaturation`))
    
    pp <- ggplot(class_data, aes(x = `Total Carbon`, 
                                 y = `Total Unsaturation`,
                                 size = -log10(`p-value`)))+
      geom_point(stat = "identity",aes(colour = `log2FoldChange`))+
      geom_point(stat = "identity",colour="black",shape=21,alpha=0.5)+
      scale_color_gradient2(low = "navy",mid="#f1f7f3",high = "firebrick3",midpoint=0,
                            limits=c(floor(min(class_data$`log2FoldChange`)),
                                     ceiling(max(class_data$`log2FoldChange`))))+
      labs(title = Class[k],colour = "Log2(FC)",size = "-log10(Pvalue)",x="Carbon",y="Degree of unsaturation")+
      theme_bw()+
      guides(size=guide_legend(order=1))+
      theme(plot.title = element_text(hjust = 0.5,size = plot.caption_size),
            plot.caption = element_text(hjust = 0.5),
            axis.title=element_text(size = axis.title_size,family=family),
            axis.text.x=element_text(size = axis.text_size,family=family),
            axis.text.y=element_text(size = axis.text_size,family=family),
            axis.title.x=element_text(size = axis.title_size,family=family,vjust=0.2),
            axis.title.y=element_text(size = axis.title_size,family=family,vjust=1),
            legend.title = element_text(size=legend.title_size,family=family),
            legend.text=element_text(size=legend.text_size,family=family),
            legend.key.size = unit(0.4, "cm"),
            plot.margin = unit(c(2,2,2,2), "cm"))
    
    ggplotsave(plot = pp,
               savepath = savepath,
               mapname = paste0(mapname,"-",Class[k]),
               width = width,
               height = height,
               dpi = dpi,
               family = family,
               ...)
    
  }
  
  if(wbsave){ 
    savexlsx1(data = result,
              filename = paste0(savepath,"/",savename))
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_lipid_point <- map_autodraw$new(moudle = Lipid_point_plot)$draw
  
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
  
  # parser$add_argument("-pd","--plotdata",default = "plotdata.xlsx", help = "原始数据")
  # parser$add_argument("-gd","--groupdata",default = "group.xlsx", help = "分组信息")
  # parser$add_argument("-wb","--wbsave",default = T, help = "是否导出绘图数据")
  # parser$add_argument("-it","--imagetype",default =c("jpg", "pdf"), help = "图片保存格式")
  # parser$add_argument("-pcs","--plot.caption_size",default = 19, help = "标题字体大小")
  # parser$add_argument("-ates","--axis.text_size",default = 14, help = "轴标签字体大小")
  # parser$add_argument("-atis","--axis.title_size",default = 16, help = "轴标题字体大小")
  # parser$add_argument("-ltis","--legend.title_size",default = 15 ,help = "图例标题字体大小")
  # parser$add_argument("-ltes","--legend.text_size",default = 10, help = "图例文本字体大小")
  # parser$add_argument("-f","--family",default = "sans", help = "文本字体")
  # parser$add_argument("-w","--width",default = 12, help = "图片保存宽度")
  # parser$add_argument("-he","--height",default = 8, help = "图片保存高度")
  # parser$add_argument("-dpi","--dpi",default = 300, help = "分辨率")
  # parser$add_argument("-sdp","--savedatapath",default = "./", help = "数据保存路径")
  # parser$add_argument("-fp","--frompath",default = "./", help = "数据来源路径")
  # parser$add_argument("-sn","--savename",default = "绘图数据-point.xlsx", help = "绘图数据保存名称")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  result <- do.call(what =  map_lipid_point,args = args) 
} 

#' @export
map_lipid_point <- map_autodraw$new(moudle = Lipid_point_plot)$draw
