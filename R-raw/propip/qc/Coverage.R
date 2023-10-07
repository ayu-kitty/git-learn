#!/opt/conda/bin/Rscript
#' 覆盖度饼图
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
seqcoverage <- function(savepath="./Statistics/", inputpath="./",inputfile="Protein quantitation.xlsx", 
                             imagetype=c("pdf","png"), height=9, width=10, dpi=300, fontfamily="sans", ...){
  pacman::p_load(ggplot2,dplyr)
  createdir(filename = savepath)
  cols<-c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#FDBF6F","#E31A1C")
  data_pro <- readdata(paste0(inputpath,inputfile))
  if(select(data_pro,ends_with("[%]"))%>%ncol()!=0){
    ##########数据整理##########
    data <- select(data_pro,ends_with("[%]"))
    ca <- as.data.frame(matrix(ncol = 2,nrow = 7))
    names(ca) <- c("Coverage","count")
    for (i in 1:6) {
      paste0("(",10*i - 10,",",10*i,"]") -> ca[i,1]
      length(which(data > 10*i - 10 & data <= 10*i)) -> ca[i,2]
    }
    ">60" -> ca[7,1]
    length(which(data >60)) -> ca[7,2]
    ca[,3]<-paste(ca[,1],paste(round(ca$count/sum(ca$count) * 100, 1), '%', sep = ''),sep=", ")
    names(ca)[3]<-"cova"
    ########绘图########
    ca[,3] <- factor(ca[,3],levels = rev(ca[,3]))
    ca_plot<-ggplot(ca, aes(x = "sub", y = count, fill = cova)) + 
      geom_bar(stat = 'identity', position = 'stack',width = 1,col="white") +
      scale_fill_manual(values=rev(cols[1:nrow(ca)]),guide = guide_legend(reverse = T))+
      coord_polar(theta = 'y')+
      labs(x = '', y = '', title = '',fill="")+
      theme(panel.grid=element_blank(),panel.background = element_rect(fill = NA),
            axis.title= element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
            legend.text=element_text(size=15),aspect.ratio=1,
            plot.margin = unit(rep(1,4),"lines"))
    #########保存数据#########
    names(ca)[1:2]<- c("range","count")
    savexlsx(ca[,-3],paste0(savepath,"Coverage_distribution.xlsx"))
    ggplotsave(plot = ca_plot,savepath=savepath,
               mapname = "Coverage_distribution",
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily,
               ...)
  }else print("~没有覆盖度参数,不绘制饼图！")
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
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Statistics/", help = "输出结果路径，默认./Statistics/")
  args <- parser$parse_args()
  coveragedist <- do.call(seqcoverage,args = args)
}

