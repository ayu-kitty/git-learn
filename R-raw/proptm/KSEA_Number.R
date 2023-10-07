#!/opt/conda/bin/Rscript
#' 激酶基因对应的底物数量排序top20
#'
#' @param savepath 绘图文件读取及生成文件保存路径
#' @param incompare 比较组名称
#' @param imagetype 图片格式
#' @param m.cutoff 底物数量
#' @param p.cutoff 激酶富集p值
#' @param height 
#' @param width 
#' @param dpi 
#' @param fontfamily 
#' @param ... 
#' @export
KSEA_Number<-function(savepath = "KSEA_Activity/",incompare = "A-vs-B", imagetype = c("pdf", "png"),
                      m.cutoff=2,p.cutoff=0.05 , height = 8, width = 14, dpi =300, fontfamily = "sans", ...){
  pacman::p_load(dplyr,ggplot2,stringr,Hmisc)
  
  path<-paste0(savepath,incompare)
  KSLinks=readxlsx(paste0(path,"/",'Kinase-Substrate Links.xlsx'))
  number_kin_sub1<- data.frame(table(KSLinks$Kinase.Gene))
  names(number_kin_sub1)<- c("Kinase","Sub_Number")
  if (nrow(number_kin_sub1)>=20) {number_kin_sub=number_kin_sub1[order(-number_kin_sub1$Sub_Number),][1:20,]}else{
    number_kin_sub=number_kin_sub1
  }
  savexlsx(number_kin_sub,paste0(path,"/KSEA_Substrate-Number-",incompare,".xlsx"))
  p<- ggplot(number_kin_sub,aes(x=reorder(Kinase,-Sub_Number),y=Sub_Number))+
    geom_bar(stat='identity',position = position_dodge(0.5),width = 0.5,fill="#4f94cd")+
    labs(y ='Number of Substrate',x='Kinases Gene',title=paste0("Substrate Number of Each Kinase (",incompare,")"))+
    geom_text(aes(label=Sub_Number), size=4,position = position_dodge(width = 0.5), vjust=-0.5)+
    theme_bw()+
    theme(text = element_text(family = "sans",colour = "black"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size=18,hjust = 0.5),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text.x = element_text(color = "black",size=12,angle = 60,hjust = 1),
          axis.text.y = element_text(color = "black",size=12,margin = margin(0,0,0,0.5,"cm")),
          plot.margin=unit(c(2,2,2,2), "lines"),
          aspect.ratio = 4/5)
  ggsave(paste0(path ,"/KSEA_Substrate-Number-",incompare,".png"),width = 10,height = 8,dpi = 300,p)
  ggsave(paste0(path ,"/KSEA_Substrate-Number-",incompare,".pdf"),width = 10,height = 8,dpi = 300,p)
  ggplotsave(plot = p,savepath=path,
             mapname = paste0("/KSEA_Substrate-Number-",incompare),
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
  parser$add_argument("-ic","--incompare",default = "A_B", help = "比较组")
  parser$add_argument("-s","--savepath",default = "KSEA_Activity/", help = "绘图文件及图片保存路径，默认KSEA_Activity/")
  parser$add_argument("-m","--m.cutoff",default=2, type= "double",help="激酶对应底物的最小数量,默认2")
  parser$add_argument("-v","--p.cutoff",default=0.05, type= "double",help="富集p值,默认0.05")
  args <- parser$parse_args()
  KSEA_Number <- do.call(KSEA_Number,args = args)
}
