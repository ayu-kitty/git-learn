#!/opt/conda/bin/Rscript
#' 肽段位点统计图
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
pepsite <- function(savepath="./Statistics/", inputpath="./",inputfile="Peptides.xlsx", 
                    imagetype=c("pdf","png"), height=10, width=10, dpi=300, fontfamily="sans", col = NULL,...){
  pacman::p_load(ggplot2,dplyr,ggpubr)
  createdir(filename = savepath)
  #输入文件####
  peps_pre <- readdata(paste0(inputpath,inputfile))
  name <- c("Deamidation 18O (N)","Acetyl (K)","GlyGly (K)","Phospho (STY)")
  names(name) <- c("Deamidated","Acetylated","Ubiquitinated","Phosphorylated")
  proty <- names(name[name%in%colnames(peps_pre)])
  peps <- peps_pre[name[name%in%colnames(peps_pre)]]
  #数据统计####
  peps <- table(peps[peps>0]) %>% as.data.frame()
  colnames(peps) <- c("Sites", "Counts")
  if(nrow(peps)>3){
    ps<-as.data.frame(matrix(ncol = 2,nrow = 4))
    names(ps)<-c("Sites","count")
    for (i in 1:3) {
      paste0(proty,"_sites ",i," (",round(peps[i,2]*100/sum(peps[,2]),digits=2),"%)")->ps[i,1]
      peps[i,2]->ps[i,2]
    }
    sum(peps[4:nrow(peps),2])->ps[4,2]
    paste0(proty,"_sites >3 (",round(ps[4,2]*100/sum(peps[,2]),digits=2),"%)")->ps[4,1]
  }else{
    ps<-as.data.frame(matrix(ncol = 2,nrow = 4))
    names(ps)<-c("Sites","count")
    for (i in 1:3) {
      paste0(proty,"_sites ",i," (",round(peps[i,2]*100/sum(peps[,2]),digits=2),"%)")->ps[i,1]
      peps[i,2]->ps[i,2]
    }
  }
  #颜色自定义####
  if(is.null(col)){
    cols<-c("#4298b5","#dd5f32", "#e19d29", "#92b06a")
  }else{
    cols <- unlist(strsplit(col,","))
  }
  ########绘图########
  ps$Sites <- factor(ps$Sites,levels = rev(ps$Sites))
  ps_plot<-ggdonutchart(ps, "count",
                        label = "Sites",                               
                        fill = "Sites",                            
                        color = "white"
  )+
    scale_fill_manual(values=rev(cols),guide = guide_legend(reverse = T))+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),legend.text=element_text(size=14))+
    theme(legend.position = c(0.5, 0.5),legend.justification = c(0.5,0.5))+
    labs(fill="")+
    theme(plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"))
    
  #########保存数据#########
  
  savexlsx(ps,paste0(savepath,proty,"_peptide_sites.xlsx"),sheet = "sheet1")
  ggplotsave(plot = ps_plot,savepath=savepath,
             mapname = paste0(proty,"_peptide_sites"),
             width = width,
             height = height,
             imagetype = imagetype,
             family=fontfamily,...)  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 10, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-c","--col",default = NULL, help = "颜色列表，以,分隔")
  parser$add_argument("-if","--inputfile", default = "Peptides.xlsx",help = "定量结果文件名，默为Peptides.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Statistics/", help = "输出结果路径，默认./Statistics/")
  args <- parser$parse_args()
  Pepsite <- do.call(pepsite,args = args)
}
