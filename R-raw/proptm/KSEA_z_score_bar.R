#!/opt/conda/bin/Rscript
#' 激酶z-score柱状图
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
KSEA_z_score_bar<-function(savepath = "KSEA_Activity/",incompare = "A-vs-B", inputdata='KSEA Kinase Scores.xlsx' , imagetype = c("pdf", "png") ,m.cutoff=2,
                      p.cutoff=0.05 , height = 10, width = 14, dpi =300, fontfamily = "sans", ...){
  pacman::p_load(dplyr,ggplot2,stringr,Hmisc,openxlsx)

  path<-paste0(savepath,incompare)
  KSScores<-read.xlsx(paste0(path,'/',inputdata))
  KSScores1<- KSScores[KSScores$m>= m.cutoff,]
  if(nrow(KSScores1)>1){
  KSScores1 = KSScores1[order(KSScores1$z.score),]
  KSScores1$color = "#060000"
  KSScores1[(KSScores1$p.value < p.cutoff) & (KSScores1$z.score < 0), 'color'] = '#36529c'
  KSScores1[(KSScores1$p.value < p.cutoff) & (KSScores1$z.score > 0), 'color'] = '#b22221'
  KSScores1$type = "no effect"
  KSScores1[(KSScores1$p.value < p.cutoff) & (KSScores1$z.score < 0), 'type'] = 'negative'
  KSScores1[(KSScores1$p.value < p.cutoff) & (KSScores1$z.score > 0), 'type'] = 'positive'

  KSScores1 <- KSScores1 %>% mutate(
    Kinase.Gene = factor(Kinase.Gene,levels = Kinase.Gene[order(z.score)]))
  if(any(KSScores1$z.score<0)){
    a1<- KSScores1$z.score[which(KSScores1$z.score<0)]
  }else{a1<- 0}
  if(any(KSScores1$z.score>0)){
    a2<- KSScores1$z.score[which(KSScores1$z.score>0)]
  }else{a2<- 0}


      tmp = with(KSScores1, labeling::extended(min(a1),max(a2), m = 5,only.loose = TRUE))
      lm = tmp[c(1,length(tmp))]
      p<- ggplot(KSScores1, aes(y=Kinase.Gene, x=z.score))+
        geom_bar(stat='identity',width=0.5,fill=KSScores1$color)+
        scale_x_continuous(limits = lm,breaks = tmp,labels = abs(tmp))+
        labs(title = incompare)+
        theme(aspect.ratio=9/16,
              panel.grid =element_blank(),#去除网格线
              panel.background = element_blank(),#背景颜色
              axis.ticks.y = element_blank(), #去除y轴
              axis.line.x = element_line(color = "black"),
              axis.text.x = element_text(size=20,color="black"),
              axis.text.y = element_text(margin = margin(t = 2, r = 2, b = 2, l = 10),
                                         color="black",size=ifelse(nrow(KSScores1)>3,ifelse(nrow(KSScores1)<7,150/nrow(KSScores1),350/nrow(KSScores1)),40)),
              axis.title = element_text(size=24,color="black"),
              plot.title = element_text(size=30,hjust=0.5,color="black"),
              plot.margin =unit(c(2,2,2,2), "lines"))

      ggplotsave(plot = p,savepath=path,
                 mapname = paste0("/KSEA-z_score-Bar"),
                 width = width,
                 height = height,
                 imagetype = imagetype,
                 family=fontfamily,
                 ...)
    }else{
      savetxt(data = "符合筛选要求的激酶数量低于2，不作z-score打分对比图",
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
  parser$add_argument("-da","--inputdata",default = "KSEA Kinase Scores.xlsx", help = "输入绘图文件，默认KSEA Kinase Scores.xlsx")
  parser$add_argument("-s","--savepath",default = "KSEA_Activity/", help = "绘图文件及图片保存路径，默认KSEA_Activity/")
  parser$add_argument("-m","--m.cutoff",default=2, type= "double",help="激酶对应底物的最小数量,默认2")
  parser$add_argument("-v","--p.cutoff",default=0.05, type= "double",help="富集p值,默认0.05")
  args <- parser$parse_args()
  KSEA_z_score_bar <- do.call(KSEA_z_score_bar,args = args)
}
