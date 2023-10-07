#!/opt/conda/bin/Rscript
#' 肽段长度分布图
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
peplengthplot <- function(savepath="./Statistics/", inputpath="./",inputfile="Peptides.xlsx", 
                          imagetype=c("pdf","png"), height=9, width=14, dpi=300, fontfamily="sans", ...){
  pacman::p_load(openxlsx,ggplot2,dplyr)
  createdir(filename = savepath)
  data_pep <- readdata(paste0(inputpath,inputfile))
  ##########数据整理##########
  if(length(intersect(c("Annotated Sequence","Sequence","StrippedSequence"),colnames(data_pep)))<1){
    print("~**********************************************************************")
    print("~请确认肽段表中肽段序列数据的列名！")
    print("~**********************************************************************")
  }else {
    data_pep_new <- data.frame(data_pep[,intersect(c("Annotated Sequence","Sequence","StrippedSequence"),colnames(data_pep))])
    colnames(data_pep_new) <- "Sequence"
    if(strsplit(data_pep_new$Sequence[1],"")[[1]][1]=="["){
      data_pep_new$Sequence<-gsub("\\[[[:alpha:]]+\\]","",data_pep_new$Sequence)%>%gsub("-","",.)%>%gsub("[.]","",.)
    }
    data_pep_new <- data.frame(Sequence=unique(data_pep_new$Sequence))
    lengpep<-nchar(data_pep_new$Sequence)%>%table%>%as.data.frame()
    names(lengpep)<-c("pepp","length")
    lengpep$pepp<-as.numeric(as.character(lengpep$pepp))
    if(nrow(lengpep)>25){
      sum(lengpep[26:nrow(lengpep),2])->la
      pepp<-lengpep[1:25,]
      paste0(">",lengpep[25,1])->pepp[26,1]
      la->pepp[26,2]
    }else lengpep->pepp 
    pepp[,1]<-factor(pepp[,1],levels = (pepp[,1])) 
    ########绘图########
    pep_plot <- ggplot(pepp,aes(pepp,length))+
      geom_bar(stat = "identity",fill="steelblue3",colour="black",width=0.8)+
      ylim(0,max(pepp$length)*1.1)+
      labs(x=" \nPeptide Length",y="Peptide Number\n ")+
      geom_text(aes(label=length), vjust=-0.5,size=5)+
      theme_bw() + 
      theme(panel.grid=element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
      theme(plot.title = element_text(hjust = 0.5),aspect.ratio=9/16,text=element_text(size=15),
            axis.text.x = element_text(size=15,color="black"),
            axis.text.y = element_text(size = 15,color="black"),
            title = element_text(size=18,color="black"),
            plot.margin=unit(c(2,2,2,2), "lines"))
    #########保存数据#########
    names(pepp)<-c("Peptides_length","Count")
    savexlsx(pepp,paste0(savepath,"Peptide_length.xlsx"))
    ggplotsave(plot = pep_plot,savepath=savepath,
               mapname = "Peptide_length",
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
  parser$add_argument("-wi","--width",default = 14, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 9, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-if","--inputfile", default = "Peptides.xlsx",help = "肽段结果文件名，默认为Peptides.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Statistics/", help = "输出结果路径，默认./Statistics/")
  args <- parser$parse_args()
  peplength <- do.call(peplengthplot,args = args)
}
