#!/opt/conda/bin/Rscript
#' 功能注释统计图
#'
#' @param savepath 结果保存路径
#' @param inputpath 输入文件路径
#' @param dbpath 背景文件路径
#' @param imagetype 保存图片类型
#' @param height 图片高度
#' @param width 图片宽度
#' @param dpi 分辨率
#' @param fontfamily 字体样式
#' @param ... 
#'
#' @export
annorateplot <- function(savepath="./Annotation/", inputpath="./",dbpath="../background/",
                         imagetype=c("pdf","png"), height=9, width=10, dpi=300, fontfamily="sans", ...){
  pacman::p_load(ggplot2,dplyr,scales)
  createdir(filename = savepath)
  profile<-dir(inputpath,pattern = "Pro.*")[1]
  sitefile<-dir(inputpath,pattern = "*Sites.*")[1]
  if(!is.na(sitefile)){
    prod<-readdata(paste0(inputpath,"/",sitefile))
  }else prod<-readdata(paste0(inputpath,"/",profile))
  pro<-data.frame(unique(prod$Accession))
  names(pro)[1]<-"Accession"
  annfile<-dir(dbpath,pattern = "Annotation.*.xlsx$")[1]
  anno<-readdata(paste0(dbpath,annfile))[,-1]
  names(anno)[1]<-"Accession"
  annopro<-merge(pro,anno,by="Accession",all.x = T)
  annopro<-annopro[,colSums(!is.na(annopro))!=0,drop = F]
  savexlsx(annopro,paste0(savepath,"annotation.xlsx"))
  annopro[annopro=="--"]<-NA
  #统计注释率
  type <- c("GO","pathway","eggNOG","PFAM","InterPro","WikiPathway","Reactome")
  names(type) <- c("GO","KEGG","eggNOG","PFAM","InterPro","WikiPathways","Reactome")
  annodata <- data.frame("Database"=names(type))
  for(i in 1:length(type)){
    floc<-grep(type[i],names(annopro))
    if(length(floc)!=0){
      annodata[i,2]<-length(na.omit(annopro[,floc[1]]))
    }
  }
  annodata[,3]<-annodata[,2]/nrow(annopro)
  names(annodata)[2:3]<-c("Number","Rate")
  annodata<-na.omit(annodata)
  if(nrow(annodata)==7){
    width=10
  }else width=ceiling(nrow(annodata)/7*10)+1
  ##########绘图##########
  colorall <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF")
  names(colorall) <-names(type)
  color <- colorall[annodata[,1]]
  annodata[,1] <- factor(annodata[,1],levels = annodata[,1])
  
  ########绘图########
  anno_plot <- ggplot(annodata,aes(Database,Rate,fill=Database))+
    geom_bar(stat = "identity",width=0.5)+
    geom_text(mapping = aes(label = annodata[,2]),vjust = -0.5,size=5)+
    labs(title = "Functional annotation of Proteins",x="", y="Percent\n ",fill="") +
    scale_fill_manual(values=color) +
    theme_bw() + 
    theme(panel.grid=element_blank())+
    theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1))+
    theme(text=element_text(size=15),axis.text.x=element_text(angle=45,hjust = 1,size=14,color="black"),
          axis.text.y=element_text(size=14,color="black"),
          title = element_text(size = 18,color="black"),
          plot.margin = unit(rep(2,4),"lines"),plot.title = element_text(hjust = 0.5))
  
  #########保存数据#########
  savexlsx(annodata,paste0(savepath,"Annotation_statistics.xlsx"))
  ggplotsave(plot = anno_plot,savepath=savepath,
             mapname = "Annotation_statistics",
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
  parser$add_argument("-wi","--width",default = 10, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 9, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-sh","--dbpath", default = "../background/",help = "背景文件存放位置,默认为../background/")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Annotation/", help = "输出结果路径，默认./Annotation/")
  args <- parser$parse_args()
  annorate <- do.call(annorateplot,args = args)
}
