#!/opt/conda/bin/Rscript
#' 肽段数分布直方图
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
pepcountplot <- function(savepath="./Statistics/", inputpath="./",inputfile="Protein quantitation.xlsx", 
                         imagetype=c("pdf","png"), height=9, width=12, dpi=300, fontfamily="sans", ...){
  pacman::p_load(openxlsx,ggplot2,dplyr)

  data_pro <- readdata(paste0(inputpath,"/",inputfile))
  ######数据处理######
  if(length(intersect(c("Peptides","# Peptides"),names(data_pro)))<1){
    print("~**********************************************************************")
    print("~蛋白表中肽段数目的列名既不是Peptides也不是# Peptides！")
    print("~**********************************************************************")
  }else {
    pepn<-select(data_pro,c(1,intersect(c("Peptides","# Peptides"),names(data_pro))))
    colnames(pepn) <- c("Accession", "Freq")
    mypep <- pepn
    if(max(mypep$Freq)>20){
      pepp<-as.data.frame(matrix(ncol = 2,nrow = 21))
      names(pepp)<-c("pepp","count")
      for (i in 1:20) {
        i->pepp[i,1]
        length(which(mypep$Freq==i))->pepp[i,2]
      }
      ">20"->pepp[21,1]
      length(which(mypep$Freq>15))->pepp[21,2]
    }else{
      pepp<-as.data.frame(matrix(ncol = 2,nrow = max(mypep$Freq)))
      names(pepp)<-c("pepp","count")
      for (i in 1:max(mypep$Freq)) {
        i->pepp[i,1]
        length(which(mypep$Freq==i))->pepp[i,2]
      }
    } 
    pepp[,1]<-factor(pepp[,1],levels = (pepp[,1])) 
    ########绘图########
    pep_plot <- ggplot(pepp,aes(pepp,count))+
      geom_bar(stat = "identity",fill="steelblue3",colour="black",width=0.8)+
      ylim(0,max(pepp$count)*1.1)+
      labs(x=" \nPeptide Number",y="Number of identified proteins\n ")+
      geom_text(aes(label=count), vjust=-0.5,size=5)+
      theme_bw() + 
      theme(panel.grid=element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
      theme(axis.text.x = element_text(size=15,color="black"),
            axis.text.y = element_text(size = 15,color="black"),
            title = element_text(size=18,color="black"),
            plot.title = element_text(hjust = 0.5),text=element_text(size=15),
            plot.margin=unit(c(2,2,2,2), "lines"))
    #########保存数据#########
    names(pepp)<- c("range","count")
    savexlsx(pepp,paste0(savepath,"/Peptide_number.xlsx"))
    ggplotsave(plot = pep_plot,savepath=savepath,
               mapname = "Peptide_number",
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily,
               ...)
    
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 12, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 9, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-if","--inputfile", default = "Protein quantitation.xlsx",help = "定量结果文件名，默认为Protein quantitation.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Statistics/", help = "输出结果路径，默认./Statistics/")
  args <- parser$parse_args()
  pepcount <- do.call(pepcountplot,args = args)
}

