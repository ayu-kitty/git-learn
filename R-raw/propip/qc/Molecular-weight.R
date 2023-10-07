#!/opt/conda/bin/Rscript
#' 蛋白分子质量分布图
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
MWplot <- function(savepath="./Statistics/", inputpath="./",inputfile="Protein quantitation.xlsx", 
                   imagetype=c("pdf","png"), height=9, width=14, dpi=300, fontfamily="sans", ...){
  
  pacman::p_load(ggplot2,dplyr)
  createdir(filename = savepath)
  data_pro <- readdata(paste0(inputpath,inputfile))
  ##########数据整理##########
  data <- select(data_pro,ends_with("[kDa]"))
  data <-na.omit(data)
  if ( max(data) > 200) {
    mw <- as.data.frame(matrix(ncol = 2,nrow = 21))
    names(mw) <- c("weight","count")
    for (i in 1:20) {
      paste0("(",10*i - 10,",",10*i,"]") -> mw[i,1]
      length(which(data > 10*i - 10 & data <= 10*i)) -> mw[i,2]
    }
    ">200" -> mw[21,1]
    length(which(data > 200)) -> mw[21,2]
  }else{
    ceiling(max(data)/10) -> num
    mw <- as.data.frame(matrix(ncol = 2,nrow = num))
    names(mw) <- c("weight","count")
    for (i in 1:num) {
      paste0("(",10*i - 10,",",10*i,"]") -> mw[i,1]
      length(which(data > 10*i - 10 & data <= 10*i)) -> mw[i,2]
    }
  }
  mw[,1] <- factor(mw[,1],levels = (mw[,1]))
  ########绘图########
  mw_plot <- ggplot(mw,aes(weight,count))+
    geom_bar(stat = "identity",fill="steelblue3",colour="black",width=0.8)+
    labs(x=" \nMolecular weight(kDa)", y="Number of identified Proteins\n ") +
    ylim(0,max(mw$count)*1.1) +
    geom_text(mapping = aes(label = paste(count)),vjust = -0.5)+
    theme_bw() + 
    theme(panel.grid=element_blank())+
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
    theme(text=element_text(size=15),#aspect.ratio=9/16,
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size=15,color="black"),
          axis.text.y=element_text(size=15,color="black"),
          title = element_text(size=18,color="black"),
          plot.margin = unit(rep(2,4),"lines"))
  
  #########保存数据#########
  names(mw)<- c("range","count")
  savexlsx(mw,paste0(savepath,"Molecular_weight.xlsx"))
  ggplotsave(plot = mw_plot,savepath=savepath,
             mapname = "Molecular_weight",
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
  parser$add_argument("-wi","--width",default = 14, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 9, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-if","--inputfile", default = "Protein quantitation.xlsx",help = "定量结果文件名，默认为Protein quantitation.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Statistics/", help = "输出结果路径，默认./Statistics/")
  args <- parser$parse_args()
  MW <- do.call(MWplot,args = args)
}

