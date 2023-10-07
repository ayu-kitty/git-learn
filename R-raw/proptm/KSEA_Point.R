#!/opt/conda/bin/Rscript
#' 激酶活性气泡图
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
KSEA_point<-function(savepath = "KSEA_Activity/",incompare = "A-vs-B", imagetype = c("pdf", "png") ,
                     m.cutoff=2,p.cutoff=0.05 , height = 10, width = 8, dpi =300, fontfamily = "sans", ...){
  pacman::p_load(dplyr,ggplot2,stringr,Hmisc)

  path<-paste0(savepath,incompare)
  KSScores<-readxlsx(paste0(path,"/","KSEA Kinase Scores.xlsx"))
  KSScores1<- KSScores[KSScores$m>= m.cutoff,]
  if(nrow(KSScores1)>0){
    if (nrow(KSScores1)>=20) {drawdata=KSScores1[order(KSScores1$p.value),][1:20,]}else{
      drawdata=KSScores1
    }
  KSLinks=readxlsx(paste0(path,"/",'Kinase-Substrate Links.xlsx'))
  da_KSlinks<- KSLinks[which(KSLinks$Kinase.Gene%in%drawdata$Kinase.Gene),]
  da_KSlinks$FC<- 2^da_KSlinks$log2FC
  da_KSlinks$change<- ifelse(da_KSlinks$FC>1,"Up","Down")
  da_KSlinks$Kinase.Gene <- factor(da_KSlinks$Kinase.Gene, levels=drawdata$Kinase.Gene)
  savexlsx(da_KSlinks,paste0(path,"/KSEA_point-",incompare,".xlsx"))
  a = length(unique(da_KSlinks$Kinase.Gene))
  p<- ggplot(da_KSlinks,aes(x=Kinase.Gene,y=paste0(`Substrate.Gene`,"_",`Substrate.Mod`),color=change))+
    geom_point(aes(size=abs(log2FC)),alpha=0.75, shape=16)+
    labs(y ='Substrates Gene',x='Kinase Gene',title= incompare,color="",size='log2(FC)')+
    scale_color_manual(values =c('#36529c', "#b22221"))+
    guides(size = guide_legend(order =1)) +
    theme_bw()+
    theme(text = element_text(family = "sans",colour = "black"),
          panel.border=element_rect(colour = "black"),
          plot.title = element_text(size=22,hjust=0.5),
          axis.text.x = element_text(size=ifelse(a>3,ifelse(a<7,80/a,250/a),20),
                                     angle = 60,hjust = 1,colour = "black"),
          axis.text.y = element_text(colour = "black",size=ifelse(nrow(da_KSlinks)>25,580/nrow(da_KSlinks),ifelse(nrow(da_KSlinks)<10,50/nrow(da_KSlinks),260/nrow(da_KSlinks)))),
          axis.title=element_text(size=15),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          legend.key.size=unit(20,'pt'),
          plot.margin=unit(c(2,2,2,2), "lines"),
          aspect.ratio = 5/3)
    ggplotsave(plot = p,savepath=path,
               mapname = paste0("/KSEA-point-",incompare),
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily,
               ...)
  }else{
    savetxt(data = "激酶底物数量低于2，不绘制气泡图",
            filename = paste0(path,"/说明.txt"),append = T)
    return()
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 8, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")

  # 此图参数
  parser$add_argument("-ic","--incompare",default = "A_B", help = "比较组")
  parser$add_argument("-s","--savepath",default = "KSEA_Activity/", help = "绘图文件及图片保存路径，默认KSEA_Activity/")
  parser$add_argument("-m","--m.cutoff",default=2, type= "double",help="激酶对应底物的最小数量,默认2")
  parser$add_argument("-v","--p.cutoff",default=0.05, type= "double",help="富集p值,默认0.05")
  args <- parser$parse_args()
  KSEA_point <- do.call(KSEA_point,args = args)
}
