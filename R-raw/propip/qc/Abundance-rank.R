#!/opt/conda/bin/Rscript
#' 蛋白表达丰度图
#'
#' @param savepath 结果保存路径
#' @param inputpath 输入文件路径
#' @param inputfile 输入文件名称
#' @param sample sample表名称
#' @param imagetype 保存图片类型
#' @param height 图片高度
#' @param width 图片宽度
#' @param dpi 分辨率
#' @param fontfamily 字体样式
#' @param ... 
#'
#' @export
abundrankplot <- function(savepath="./Statistics/", inputpath="./",inputfile="Protein quantitation.xlsx",
                          sample = "Sample_Group.xlsx",imagetype=c("pdf","png"), height=9, width=10, dpi=300, fontfamily="sans", ...){

  pacman::p_load(ggplot2,dplyr,tidyr,ggrepel)
  createdir(filename = savepath)
  data_pro <- readdata(paste0(inputpath,inputfile))
  sam1 <- readdata(paste0(inputpath,"/",sample))
  sam1 <- sam1[!duplicated(sam1[,1]),]
  rownames(sam1)<-sam1[,1]
  thisam<-intersect(sam1[,1],names(data_pro))
  sam<-sam1[thisam,]
  if("ID" %in% colnames(data_pro)){
    data_pro["Label"]<-data_pro["ID"]
  }else data_pro["Label"]<-data_pro["Accession"]
  data <- select(data_pro,c("Label",all_of(sam[,1])))
  
  data$Median_abundance <- apply(as.data.frame(data[,-1]),1,median,na.rm = T)
  dataquantity <- data[order(-data$Median_abundance),]
  savexlsx(dataquantity,paste0(savepath,"Abundances_rank.xlsx"))#输出表格文件

  if(nrow(dataquantity)>10){
    datapre <- dataquantity %>% drop_na(Median_abundance)
    dataplot <- data.frame(name=datapre$Label,number=1:nrow(datapre),quantity=datapre$Median_abundance)
    dataplot$sort <- ifelse(dataplot$number %in% 1:5,"top","normal")
    dataplot<-arrange(dataplot,sort)
    dataplot$label <- ifelse(dataplot$number %in% 1:5,as.character(dataplot$name),"")
    dataplot$name <- factor(dataplot$name,levels = dataplot$name)
    ########绘图########
    quantity_plot <- ggplot(aes(x=number,y=log10(quantity),colour=sort),data = dataplot)+
      geom_point(size = 2.5)+ylim(NA,max(log10(dataplot$quantity)*1.05))+
      scale_color_manual(values = c("top" = "#c60008","normal" = "#9A9FB6"))+
      labs(x=" \nProtein rank",y="Median abundance (log10)\n ")+
      theme_bw() + 
      theme(panel.grid=element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
      geom_text_repel(aes(number, log10(quantity), label = label),color="black", 
                      size=3, segment.size=0.5, nudge_x=nrow(dataplot)*0.2, direction="y", hjust=0) +
      theme(text=element_text(size=15),aspect.ratio=9/10,legend.position = 'none',
            axis.text.x = element_text(vjust = -1,size=15,color="black"),
            axis.text.y=element_text(vjust = 0.5,size=15,color="black"),
            title = element_text(size=18,color="black"),
            plot.margin = unit(rep(2,4),"lines"))
    ggplotsave(plot = quantity_plot,savepath=savepath,
               mapname ="Abundances_rank",
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily,
               ...)
  }else{
    savetxt(data = "蛋白数量太少，不绘制表达丰度图",
            filename = paste0(savepath,"说明.txt"),append = T)
    return()
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 10, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 9, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-if","--inputfile", default = "Protein quantitation.xlsx",help = "定量结果文件名，默认为Protein quantitation.xlsx")
  parser$add_argument("-samp","--sample", default = "Sample_Group.xlsx",help = "包含样本信息的数据文件名，默认为sample_information.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Statistics/", help = "输出结果路径，默认./Statistics/")
  args <- parser$parse_args()
  abundrank <- do.call(abundrankplot,args = args)
}

