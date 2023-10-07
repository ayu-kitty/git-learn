#!/opt/conda/bin/Rscript

#' 蛋白实验标准曲线图
#'
#' @param inputpath 输入文件路径
#' @param imagetype 保存图片类型
#' @param height 图片高度
#' @param width 图片宽度
#' @param dpi 分辨率
#' @param fontfamily 字体样式
#' @param ... 
#' @export
bqplot <- function(inputpath="./rawdata/pic/",
                   imagetype=c("png"), height=6, width=10, dpi=300, fontfamily="sans",
                   ...){
  pacman::p_load(ggplot2,dplyr,ggpmisc)
  inputfile=dir(inputpath,pattern="^table.*.xlsx")[1]
  bqname<-getsheetname(paste0(inputpath,"/",inputfile))
  if("standard_concentration" %in% bqname){
    dat<-readxlsx(paste0(inputpath,"/",inputfile),sheet = "standard_concentration") %>% t() %>% as.data.frame.array()
    names(dat)<-dat[1,]
    dat<-dat[-1,]
    dat<-apply(dat, 2, as.numeric) %>% as.data.frame()
    for(i in 1:(ncol(dat)-1)){
      bq<-select(dat,c(i,ncol(dat)))
      #绘图
      p<-ggplot(bq, aes(x = bq[,1], y = bq[,2])) +
        geom_point(shape=18,color="royalblue",size=7,fill="royalblue") +
        geom_smooth(method = "lm",fill=NA,col="grey46")+
        stat_poly_eq(
          formula = y ~ x,
          use_label(c("eq")),
          parse = TRUE,size=7,label.y = 0.95
        )+
        stat_poly_eq(
          formula = y ~ x,
          use_label(c("R2")),rr.digits = 3,
          parse = TRUE,size=7,label.y = 0.85
        )+
        theme_bw() + 
        theme(panel.grid=element_blank())+
        theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1))+
        scale_x_continuous(limits=c(min(bq[,1]),ceiling(max(bq[,1])*10)/10),breaks=seq(min(bq[,1]),ceiling(max(bq[,1])*10)/10,0.1))+
        scale_y_continuous(limits=c(min(bq[,2]),ceiling(max(bq[,2])*10)/10),breaks=seq(min(bq[,2]),ceiling(max(bq[,2])*10)/10,0.1))+
        labs(x=expression("OD"[562]),y=expression(paste("Protein Concentration (",mu,"g/",mu,"L)")),title = "Protein quantitation standard curve")+
        theme(legend.text=element_text(size=15))+
        theme(plot.title = element_text(hjust = 0.5, size=20))+
        theme(axis.text.x=element_text(size=15,color="black"),axis.text.y=element_text(size=15,color="black"),title = element_text(size=18,color="black"))+
        theme(plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"))
      
      
      ggplotsave(plot = p,savepath=inputpath,
                 mapname = paste0("1.standard_curve",i),
                 width = width,
                 height = height,
                 imagetype = imagetype,
                 family=fontfamily,
                 ...)
      
    }
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 10, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 6, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  parser$add_argument("-ip","--inputpath", default = "./rawdata/pic/",help = "实验pic文件夹路径，默认为./rawdata/pic/")
  args <- parser$parse_args()
  bqplot <- do.call(bqplot,args = args)
}

