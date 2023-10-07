#!/opt/conda/bin/Rscript

#' 蛋白表达密度图
#'
#' @param savepath 保存路径
#' @param inputpath 输入文件路径
#' @param inputfile 输入文件名称
#' @param classfile classtype.xlsx
#' @param imagetype 保存图片类型
#' @param height 高度
#' @param width 宽度
#' @param dpi 分辨率
#' @param fontfamily 字体类型 
#' @param ... 
#' @export
prodensityplot <- function(savepath="./Expression/Densityplot/", inputpath="./",inputfile="绘图数据.xlsx",
                           imagetype=c("pdf","png"), height=10, width=NULL, dpi=300, fontfamily="sans", 
                           classfile = "../classtype.xlsx",
                           ...){
  pacman::p_load(ggplot2,dplyr,reshape2)
  
  dat<-readxlsx(paste0(inputpath,"/",inputfile))[,-1]
  samplecol<-stylefun_group(sample = T,classfile = classfile)[names(dat)]
  if(is.null(width)){
    width = ifelse(ncol(dat)>80,6+(ceiling(ncol(dat)/20)),9+3*(ceiling(ncol(dat)/16)))
  }
  #绘图
  plotdat<-na.omit(melt(dat))
  plotdat$variable<-factor(plotdat$variable,levels=unique(plotdat[,1]))
  p<-ggplot(plotdat,aes(x=value,color=variable))+
    stat_density(position = "identity",geom = 'line',alpha=0.5,size=1)+
    theme_bw() + 
    theme(panel.grid=element_blank())+
    theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1))+
    labs(x="",y="Density",color="") +
    scale_color_manual(values=samplecol)+
    guides(color=ifelse(length(unique(plotdat$variable))>80,"none","legend"))+
    theme(legend.text=element_text(size=15))+
    theme(axis.text.x=element_text(angle=60,hjust = 1,size=15,color="black"),axis.text.y=element_text(size=15,color="black"),title = element_text(size=18,color="black"))+
    theme(plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"))
  
  ggplotsave(plot = p,savepath=savepath,
             mapname = "Densityplot",
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
  parser$add_argument("-i","--imagetype",default = c("png","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf"))
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "图片高度")
  parser$add_argument("-cf","--classfile",default = "../classtype.xlsx", help = "模板")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-if","--inputfile", default = "绘图数据.xlsx",help = "定量结果文件名，默认为绘图数据.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Expression/Densityplot/", help = "输出结果路径，默认./Expression/Densityplot/")
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  prodensity <- do.call(prodensityplot,args = args)
}
