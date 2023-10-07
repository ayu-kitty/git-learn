#!/opt/conda/bin/Rscript
#' 样本蛋白数鉴定分布图
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
sampleproplot <- function(savepath="./Identified/", inputpath="./",inputfile="Protein quantitation.xlsx",
                          sample = "Sample_Group.xlsx",imagetype=c("pdf","png"), 
                          height=NULL, width=NULL, dpi=300, fontfamily="sans", 
                          classfile = "../classtype.xlsx",
                          ...){
  pacman::p_load(ggplot2,dplyr)
  createdir(filename = savepath)
  ##########数据整理##########
  data_pro <- readdata(paste0(inputpath,inputfile))
  if(length(unique(data_pro[,1]))!=nrow(data_pro)){
    protype="Site"
  }else protype="Protein"
  sample <- readdata(paste0(inputpath,sample),sheet = 1)
  sample <- sample[which(sample$Group != "删除"),]
  sample<-sample[!(duplicated(sample[,1])),]
  data<-data_pro[,sample[,1]]
  rownames(sample)<-sample[,1]
  data[data == 0] <- NA
  data_sample <- data.frame(sample=colnames(data),number=apply(data,2,function(x) length(na.omit(x))),group=sample[colnames(data),2])
  data_sample$group<- as.character(data_sample$group)
  data_sample$sample=factor(data_sample$sample,level=data_sample$sample)
  if(length(sample[,1])>30){
    data_sample$rank <- 1:nrow(data_sample)
    ## 多样本绘图函数
    sample_plot <- ggplot(data_sample,mapping = aes(x = rank,y = number,group = 1,color = group))+geom_line() + geom_point(size = 4) +
      scale_color_manual(values = stylefun_group(classfile = classfile)[data_sample$group])+
      labs(title = paste0("Median Number = ",round(median(data_sample$number))),x=" \nSample Rank", y=paste0("Number of ",protype,"(Sample)\n "),fill="") +
      geom_hline(aes(yintercept=median(number)), color="gray50", linetype="dashed", size=0.75)+
      theme(panel.grid=element_blank(),panel.background = element_rect(fill = NA),text=element_text(size=15),
            panel.border = element_rect(fill=NA,color="black",size = 1,linetype="solid"),plot.margin = unit(rep(2,4),"lines"),
            axis.text.x=element_text(vjust = -1,size=15,color="black"),
            axis.text.y=element_text(vjust = 0.5,size=15,color="black"),
            legend.text=element_text(size=15),
            title = element_text(size = 16,color="black"))
    height <- ifelse(is.null(height),9,height)
    width <- ifelse(is.null(width),15,width)
    
  }else{
    ## 少样本绘图函数
    sample_plot <- ggplot(data_sample,aes(sample,number,fill=group))+
      scale_fill_manual(values=stylefun_group(classfile = classfile)[data_sample$group]) +
      coord_flip()+geom_bar(stat = "identity",width=0.6)+
      labs(x="", y=paste0("\nNumber of ",protype,"(Sample)"),fill="") +
      ylim(0,max(as.numeric(data_sample[,2]))*1.2) +
      geom_text(mapping = aes(label = number),hjust = -0.5,size=5)+
      theme(panel.grid=element_blank(),panel.background = element_rect(fill = NA),text=element_text(size=15),
            panel.border = element_rect(fill=NA,color="black",linewidth = 1,linetype="solid"),plot.margin = unit(rep(2,4),"lines"),
            legend.text=element_text(size=15),
            axis.text.x=element_text(size=15,vjust = -1,color="black"),
            axis.text.y=element_text(size=15,color="black"),
            title = element_text(size = 18,color="black"))
    char <- nchar(sample$Sample)
    wi <- ifelse(max(char)>20,12,10)
    height <- ifelse(is.null(height),10,height)
    width <- ifelse(is.null(width),wi,width)
  }
  #########保存数据#########
  savexlsx(data_sample,paste0(savepath,"Sample_protein_Number.xlsx"))
  ggplotsave(plot = sample_plot,savepath=savepath,
             mapname = "Sample_protein_Number",
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
  parser$add_argument("-wi","--width",default = NULL, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = NULL, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  parser$add_argument("-cf","--classfile",default = "../classtype.xlsx", help = "模板")
  
  # 此图参数
  parser$add_argument("-if","--inputfile", default = "Protein quantitation.xlsx",help = "定量结果文件名，默认为Protein quantitation.xlsx")
  parser$add_argument("-samp","--sample", default = "Sample_Group.xlsx",help = "包含样本信息的数据文件名，默认为Sample_Group.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Identified/", help = "输出结果路径，默认./Identified/")
  args <- parser$parse_args()
  samplepro <- do.call(sampleproplot,args = args)
}
