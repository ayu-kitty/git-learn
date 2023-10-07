#!/opt/conda/bin/Rscript
#' 基本鉴定结果统计图
#'
#' @param savepath 结果保存路径
#' @param inputpath 输入文件路径
#' @param inputfile 输入文件名称
#' @param imagetype 保存图片类型
#' @param height 图片高度
#' @param width 图片宽度
#' @param dpi 分辨率
#' @param fontfamily 字体样式
#' @param ... 
#'
#' @export
identifiedplot <- function(savepath="./Statistics/", inputpath="./",inputfile="Identified_number.xlsx", mapname = "Identified-number",
                           imagetype=c("pdf","png"), label=NULL,height=10, width=10, dpi=300, fontfamily="sans", ...){
  pacman::p_load(openxlsx,ggplot2,dplyr)
  createdir(filename = savepath)
  pp <- readdata(paste0(inputpath,inputfile))
  ##########数据整理##########
  if(ncol(pp)==5){
    if(is.null(label)){
      pp<-pp[1,-1]
    }else{
      pp<-pp[pp[,1]==label,-1]
    }
    names(pp)<-c("Total spectrum","Matched spectrum","Identified proteins","Identified peptides")
    ppc<-t(pp)%>%as.data.frame()
    ppc[,2]<-rownames(ppc)
    names(ppc)<-c("count","type")
    ppc[,2]=factor(ppc[,2],levels=ppc[,2])
    ########绘图########
    pp_plot<- ggplot(ppc,aes(ppc[,2],ppc[,1],fill=ppc[,2]))+
      geom_bar(stat = "identity",width=0.5)+
      labs(x="", y="",fill="") +
      ylim(0,max(ppc[,1])*1.1) +
      theme_bw() + 
      theme(panel.grid=element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
      scale_fill_manual(values=c("#F5CB5C","#8bc24c","#1687a7","#eb5c2f")) +
      geom_text(mapping = aes(label = count),vjust = -0.5,size=5)+
      theme(axis.text.x=element_text(angle=45,hjust = 1,size=15,color="black"),axis.text.y=element_text(size=15,color="black"),
            plot.margin = unit(rep(2,4),"lines"),text=element_text(size=15),aspect.ratio=1)

  }else if(ncol(pp)==4){
    #数据统计####
    p3_data <- data.frame(V1=c("phos_pro","phos_pep","phos_site"),
                          V2=pp[1,2:4] %>% as.numeric())
    N <- p3_data$V2
    ymax <- max(N)
    names(p3_data)<-c("type","count")
    p3_data[,1]<-c("Proteins", "Peptides", "Sites")
    p3_data[,1]<-factor(p3_data[,1],levels = p3_data[,1])
    ########绘图########
    pp_plot <- ggplot(p3_data,aes(type,count,fill=type))+
      geom_bar(stat = "identity",width=0.5)+
      ylim(0,ymax*1.1)+
      labs(x="",y="",fill="")+
      theme_bw() + 
      theme(panel.grid=element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
      scale_fill_manual(values=c("Proteins"="#1687a7","Peptides"="#eb5c2f","Sites"="#41924B"))+
      geom_text(aes(label=count), vjust=-0.5,size=5)+
      theme(plot.margin=unit(c(2,2,2,2), "lines"))+
      theme(axis.text.x = element_text(size=15,color="black"),axis.text.y=element_text(size=15,color="black"),legend.text=element_text(size=15))
  }else if(ncol(pp)==3 & nrow(pp)==1){
    #数据统计####
    ppdata <- data.frame(V1=c("Proteins","Peptides"),
                          V2=pp[1,2:3] %>% as.numeric())
    ymax <- max(ppdata$V2)
    names(ppdata)<-c("type","count")
    ppdata[,1]<-factor(ppdata[,1],levels = ppdata[,1])
    ########绘图########
    pp_plot <- ggplot(ppdata,aes(type,count,fill=type))+
      geom_bar(stat = "identity",width=0.4)+
      ylim(0,ymax*1.1)+
      labs(x="",y="",fill="")+
      theme_bw() + 
      theme(panel.grid=element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
      scale_fill_manual(values=c("Proteins"="#1687a7","Peptides"="#eb5c2f"))+
      guides(color="none")+
      geom_text(aes(label=count), vjust=-0.5,size=5)+
      theme(plot.margin=unit(c(2,2,2,2), "lines"))+
      theme(axis.text.x = element_text(size=15,color="black"),axis.text.y=element_text(size=15,color="black"),legend.text=element_text(size=15))
  }else if(ncol(pp)==3 & nrow(pp)==2){
    ppdata <- data.frame(V1=c("Proteins-DIA","Peptides-DIA","Proteins-DDA","Peptides-DDA"),
                         V2=c(pp[1,2:3],pp[2,2:3]) %>% as.numeric())
    ymax <- max(ppdata$V2)
    names(ppdata)<-c("type","count")
    ppdata[,1]<-factor(ppdata[,1],levels = ppdata[,1])
    ########绘图########
    pp_plot<- ggplot(ppdata,aes(type,count,fill=type))+
      geom_bar(stat = "identity",width=0.5)+
      labs(x="", y="",fill="") +
      ylim(0,max(ppdata[,2])*1.1) +
      theme_bw() + 
      theme(panel.grid=element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
      scale_fill_manual(values=c("#F5CB5C","#8bc24c","#1687a7","#eb5c2f")) +
      geom_text(mapping = aes(label = count),vjust = -0.5,size=5)+
      theme(axis.text.x=element_text(angle=45,hjust = 1,size=15,color="black"),axis.text.y=element_text(size=15,color="black"),
            plot.margin = unit(rep(2,4),"lines"),text=element_text(size=15),aspect.ratio=1)
  }
  #########保存图片#########
  ggplotsave(plot = pp_plot,savepath=savepath,
             mapname = mapname,
             width = width,
             height = height,
             imagetype = imagetype,
             family=fontfamily,
             ...)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 10, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-if","--inputfile", default = "Identified_number.xlsx",help = "基本定量结果文件名，默认为Identified_number.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为rawdata路径")
  parser$add_argument("-mn","--mapname", default = "Identified-number",help = "图片名称，默认为Identified-number")
  parser$add_argument("-sp","--savepath",default = "./Statistics/", help = "输出结果路径，默认./Statistics/")
  args <- parser$parse_args()
  identified <- do.call(identifiedplot,args = args)
}

