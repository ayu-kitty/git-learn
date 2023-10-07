#!/opt/conda/bin/Rscript
#' 激酶活性排序top20柱状图
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
KSEA_Bar<-function(savepath = "./KSEA_Activity/",incompare = "A-vs-B", imagetype = c("pdf", "png") ,
                   m.cutoff=2,p.cutoff=0.05 , height = 8, width = 14, dpi =300, fontfamily = "sans", ...){
  pacman::p_load(dplyr,ggplot2,stringr,Hmisc)

  path<-paste0(savepath,incompare)
  KSScores<-readxlsx(paste0(path,"/","KSEA Kinase Scores.xlsx"))
  KSScores1<- KSScores[KSScores$m>= m.cutoff,]
  if(nrow(KSScores1)>1){
    if (nrow(KSScores1)>=20) {drawdata=KSScores1[order(KSScores1$p.value),][1:20,]}else{
      drawdata=KSScores1
    }
    drawdata<-drawdata[order(drawdata$z.score),]
    drawdata$Kinase.Gene<-factor(drawdata$Kinase.Gene,level=drawdata$Kinase.Gene)
    savexlsx(drawdata,paste0(path,"/KSEA_bar-",incompare,".xlsx"))
    drawdata$color = "white"
    drawdata[(drawdata$z.score < 0), ncol(drawdata)] = '#36529c'
    drawdata[(drawdata$z.score > 0), ncol(drawdata)] = '#b22221'
    p<- ggplot(drawdata,aes(x=Kinase.Gene,y=z.score,fill=z.score))+
      geom_bar(stat="identity",position = position_dodge(0.6),width = 0.6)+
      labs(y ='z.score',x='Kinase Gene',title=incompare,color="")+
      theme_bw()+
      theme(text = element_text(family = fontfamily,colour = "black"))+
      theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
            plot.title = element_text(size=18,hjust=0.5),
            axis.text.x = element_text(colour="black",hjust =1,angle = 45,size = 18),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(margin = margin(0,1,0,1,"cm"),size = 20),
            axis.ticks.y = element_line(color = "black",size = 1),
            axis.text.y = element_text(colour="black",size = 16),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 16),
            plot.margin=unit(c(2,2,2,2), "lines"),
            aspect.ratio = 2/3)
    if (unique( drawdata$z.score)%>%length(.)>1){
      p=p+scale_fill_gradient2(high='#b22221',mid = '#e0dbcd',low = '#36529c')
    }
    ggplotsave(plot = p,savepath=path,
               mapname = paste0("/KSEA-Bar-",incompare),
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily,
               ...)
  }else{
    savetxt(data = "符合筛选要求的激酶数量低于2，不作柱状图",
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
  parser$add_argument("-wi","--width",default = 14, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 9, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")

  # 此图参数
  parser$add_argument("-ic","--incompare",default = "A_B", help = "比较组")
  parser$add_argument("-s","--savepath",default = "KSEA_Activity/", help = "绘图文件及图片保存路径，默认KSEA_Activity/")
  parser$add_argument("-m","--m.cutoff",default=2, type= "double",help="激酶对应底物的最小数量,默认2")
  parser$add_argument("-v","--p.cutoff",default=0.05, type= "double",help="富集p值,默认0.05")
  args <- parser$parse_args()
  KSEA_Bar <- do.call(KSEA_Bar,args = args)
}
