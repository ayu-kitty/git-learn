#!/opt/conda/bin/Rscript
#' STY位点分布图
#'
#' @param savepath 结果保存路径
#' @param inputpath 输入文件路径
#' @param inputfile 输入文件名称
#' @param imagetype 保存图片类型
#' @param height 图片高度
#' @param width 图片宽度
#' @param dpi 分辨率
#' @param fontfamily 字体样式
#' @param col 颜色
#'
#' @export
#'
phostyplot <- function(savepath="./Statistics/", inputpath="./",inputfile="Phospho(STY) Sites.xlsx", 
                imagetype=c("pdf","png"), height=10, width=14, dpi=300, fontfamily="sans", col = NULL){
  pacman::p_load(ggplot2,dplyr,ggpubr)
  sty_pre <- readdata(paste0(inputpath,inputfile))
  sty <- table(sty_pre$`Amino acid`) %>% as.data.frame()
  colnames(sty) <- c("Site", "count")
  
  sty[1,3] <- paste("S; ", paste("Percentage:", sprintf("%.2f",sty[1,2]*100/sum(sty[,2])), "%"))
  
  sty[2,3] <- paste("T; ", paste("Percentage:", sprintf("%.2f",sty[2,2]*100/sum(sty[,2])), "%"))
  
  sty[3,3] <- paste("Y; ", paste("Percentage:", sprintf("%.2f",sty[3,2]*100/sum(sty[,2])), "%"))
  
  names(sty)[3]<-"label"
  if(is.null(col)){
    cols<-c("#EB6E44","#D3E397", "#FCE38A")
  }else{
    cols <- unlist(strsplit(col,","))
  }
  ########绘图########
  sty$label <- factor(sty$label,levels = rev(sty$label))
  sty_plot<-ggdonutchart(sty, "count",
                        label = "label",                               
                        fill = "label",                            
                        color = "white"
  )+
    scale_fill_manual(values=rev(cols),guide = guide_legend(reverse = T))+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),legend.text=element_text(size=14))+
    theme(legend.position = c(0.5, 0.5),legend.justification = c(0.5,0.5))+
    labs(fill="")+
    theme(plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"))
  #########保存数据#########
  savexlsx(sty[,2:3],paste0(savepath,"STY_locus_distribution.xlsx"),sheet = "sheet1")
  ggplotsave(plot = sty_plot,savepath=savepath,
             mapname = "STY_locus_distribution",
             width = width,
             height = height,
             imagetype = imagetype,
             family=fontfamily)  
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 14, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-c","--col",default = NULL, help = "颜色列表，以,分隔")
  parser$add_argument("-if","--inputfile", default = "Phospho (STY)Sites.xlsx",help = "定量结果文件名，默为Phospho(STY)Sites.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Statistics/", help = "输出结果路径，默认./Statistics/")
  
  args <- parser$parse_args()
  Phostyplot <- do.call(phostyplot,args = args)
}
