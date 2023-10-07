#!/opt/conda/bin/Rscript

#' 蛋白表达箱线图
#'
#' @param savepath 保存文件夹路径
#' @param inputpath 输入数据所在路径
#' @param inputfile 输入数据名称
#' @param imagetype 保存图片类型
#' @param height 高度
#' @param width 宽度
#' @param dpi 分辨率
#' @param fontfamily 字体类型
#' @param classfile classtype.xlsx
#' @param ... 
#' @export
proboxplot <- function(savepath="./Expression/Boxplot/", inputpath="./",inputfile="绘图数据.xlsx",
                       imagetype=c("pdf","png"), height=10, width=NULL, dpi=300, fontfamily="sans", 
                       classfile = "../classtype.xlsx",
                       log = F,
                       ...){
  pacman::p_load(ggplot2,dplyr,reshape2)
  dat <- readdata(paste0(inputpath,"/",inputfile))[,-1,drop = F]
  
  if(any(is.na(dat))){
    dat[is.na(dat)] <- min(as.matrix(dat),na.rm = T)
  }
  
  # if(dim(dat)[2] > 60){
  #   dat <- dat[,1:60]
  # }
  
  if(log){
    if(all(dat > 0)){
      dat <- log(dat)
    }
  }
  
  classdata <- readdata(paste0(dirname(classfile), "/classfile.yaml"))
  groupdata <- data.frame(sample = colnames(dat))
  for ( i in 1:dim(groupdata)[1]) {
    groupdata[i,"group"] <- names(classdata)[unlist(lapply(classdata, function(x,y){any(y %in% x)},y = groupdata[i,"sample"]))][1]
  }
  samplecol <- stylefun_group(classfile = classfile,subset = unique(groupdata$group))
  
  if(is.null(width)){
    width = ifelse(ncol(dat)>80,6+2*(ceiling(ncol(dat)/20)),6+6*(ceiling(ncol(dat)/16)))
  }
  #绘图
  plotdat<-na.omit(melt(dat,id.vars = 0))
  plotdat$variable <- factor(plotdat$variable,levels=unique(plotdat[,1]))
  plotdat <- merge(plotdat,groupdata,by.x = "variable",by.y = "sample")
  if(length(unique(plotdat$variable))<80){
    p<-ggplot(plotdat,aes(x=variable,y=value,color=group))+
      stat_boxplot(geom = "errorbar",width=0.35)+
      geom_boxplot(aes(fill=group))+
      geom_boxplot(alpha=0.5,outlier.size=1)+
      theme_bw() + 
      theme(panel.grid=element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
      labs(x="",y="Expression abundance (Normalization)",color="") +
      scale_color_manual(values=samplecol,breaks =names(samplecol),labels=names(samplecol))+
      scale_fill_manual(values=samplecol)+
      guides(fill="none")+
      guides(color="legend")+
      theme(axis.text.x=element_text(angle=60,hjust = 1,size=15,color="black"),
            axis.text.y=element_text(size=15,color="black"),
            title = element_text(size=18,color="black"))+
      theme(plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"))
  }else{
    p<-ggplot(plotdat,aes(x=variable,y=value,color=group))+
      stat_boxplot(geom = "errorbar",width=0.35)+
      geom_boxplot(aes(fill=group))+
      geom_boxplot(alpha=0.5,outlier.size=1)+
      theme_bw() + 
      theme(panel.grid=element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
      labs(x="",y="Expression abundance (Normalization)",color="") +
      scale_color_manual(values=samplecol,breaks =names(samplecol),labels=names(samplecol))+
      scale_fill_manual(values=samplecol)+
      guides(fill="none")+
      guides(color="legend")+
      theme(legend.text=element_text(size=15))+
      theme(axis.text.x=element_blank())+
      theme(axis.text.y=element_text(size=15,color="black"),
            title = element_text(size=18,color="black"),
            axis.ticks.length.x = unit(-0.1,"cm"))+
      theme(plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"))
  }
  ggplotsave(plot = p,savepath=savepath,
             mapname = "Boxplot",
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
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  parser$add_argument("-cf","--classfile",default = "../classtype.xlsx", help = "模板")
  parser$add_argument("-i","--imagetype",default = c("png","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf"))
  
  # 此图参数
  parser$add_argument("-if","--inputfile", default = "绘图数据.xlsx",help = "定量结果文件名，默认为绘图数据.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Expression/Boxplot/", help = "输出结果路径，默认./Expression/Boxplot/")
  parser$add_argument("-l","--log", default = F, action = "store_true", help = "是否log处理")
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  probox <- do.call(proboxplot,args = args)
}

