#!/opt/conda/bin/Rscript

#' 样本层次聚类树状图
#'
#' @param savepath 保存路径
#' @param inputpath 输入文件路径
#' @param inputfile 输入文件名
#' @param classfile classtype.xlsx
#' @param height 高度
#' @param width 宽度
#' @param fontfamily 字体
#' @param ... 
#' @export
procorrtree <- function(savepath="./Expression/Sampletreeplot/", inputpath="./",inputfile="绘图数据.xlsx",
                        height=10, width=NULL, fontfamily="sans",scaledat="F",
                        classfile = "../classtype.xlsx",
                        imagetype = c("png", "pdf"),
                        ...){
  pacman::p_load(ggplot2,dplyr,dendextend)
  dat <- readdata(paste0(inputpath,"/",inputfile))[,-1,drop = F]
  if(ncol(dat)>1){
    if(is.null(width)){
      width = ifelse(ncol(dat)>80,8+3.5*(ceiling(ncol(dat)/20)),12*ceiling(ncol(dat)/17))
    }
    #绘图
    if(is.na(sum(dat))){
      dat[is.na(dat)]<-apply(dat,2,function(x){min(na.omit(x))})%>%min()
    }
    if(scaledat=="T"){
      dat_scale <- t(apply(dat, 1, scale)) %>% as.data.frame() #数据标准化
      names(dat_scale)<-names(dat)
    }else dat_scale<-dat
    
    dam <- as.matrix(dist(t(dat_scale), method =  "euclidean"))
    savexlsx3(as.data.frame.array(dam),paste0(savepath,"/Samples_Distance.xlsx"),sheet = "Samples_Distance")
    dend <- as.dendrogram(hclust(dist(t(dat_scale), method =  "euclidean")), method = "complete")
    a <- hclust(dist(t(dat_scale), method =  "euclidean"))
    b <- a$labels[a$order]
    labels_colors(dend) <- stylefun_group(sample = T,classfile = classfile)[b]
    plotfile(savepath = savepath,mapname = "SamplesTreeplot",imagetype = imagetype, 
             width =width, height = height,family=fontfamily)
    par(mar=c(10,4,2,4) + 0.1,cex=1.5)
    plot(dend, main = NULL, sub="", xlab="")
    plotsave()
  }
  
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "图片高度")
  parser$add_argument("-cf","--classfile",default = "../classtype.xlsx", help = "模板")
  parser$add_argument("-i","--imagetype",default = c("png","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf"))
  
  # 此图参数
  parser$add_argument("-if","--inputfile", default = "绘图数据.xlsx",help = "定量结果文件名，默认为绘图数据.xlsx")
  parser$add_argument("-sl","--scaledat", default = "F",help = "数据是否经过scale，默认F，代谢需选T")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Expression/Sampletreeplot/", help = "输出结果路径，默认./Expression/Sampletree/")
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  Procorrtree <- do.call(procorrtree,args = args)
}
