#!/opt/conda/bin/Rscript
#' 蛋白位点数量分布图
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
#'
prosite <- function(savepath="./Statistics/", inputpath="./",inputfile=NULL, 
                    imagetype=c("pdf","png"), height=7.5, width=12, dpi=300, fontfamily="sans", ...){
  pacman::p_load(ggplot2,dplyr,RColorBrewer,stringr)
  createdir(filename = savepath)
  #输入文件####
  if(is.null(inputfile)){
    site <- readdata(list.files(inputpath,"*Sites.xlsx$"))
  }else{
    site <- readdata(paste0(inputpath,inputfile))
  }
  name <- paste0(c("Deamidation 18O (N)","Acetyl (K)","GlyGly (K)","Phospho (STY)")," Probabilities")
  names(name) <- c("Deamidated","Acetylated","Ubiquitinated","Phosphorylated")
  proty <- names(name[name%in%colnames(site)])  
  #数据统计####
  pro_sites_n <- table(site$Accession) %>%table()%>% as.data.frame()
  colnames(pro_sites_n) <- c("Sites", "Freq")
  if(nrow(pro_sites_n)>15){
    prop<-as.data.frame(matrix(ncol = 2,nrow = 16))
    names(prop)<-c("Sites","Count")
    for (i in 1:15) {
      pro_sites_n[i,1]->prop[i,1]
      pro_sites_n[i,2]->prop[i,2]
    }
    paste0(">",prop[15,1])->prop[16,1]
    sum(pro_sites_n[16:nrow(pro_sites_n),2])->prop[16,2]
  }else{
    prop<-pro_sites_n
    names(prop)<-c("Sites","Count")
  } 
  prop <- na.omit(prop)
  prop[,1]<-factor(prop[,1],levels = (prop[,1]))
  #绘图####
  pro_plot <- ggplot(prop,aes(Sites,Count))+
    geom_bar(stat = "identity",fill="steelblue3",colour="black",width=0.8)+
    theme_bw() + 
    theme(panel.grid=element_blank())+
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
    ylim(0,max(prop$Count)*1.1)+
    labs(x=paste0("Number of ",proty," Sites in a protein\n "),y=paste0(" \n",proty," Protein numbers"))+
    geom_text(aes(label=Count), vjust=-0.5,size=5)+
    theme(axis.text.x = element_text(size=15,color="black"),axis.text.y = element_text(size = 15,color="black"),title = element_text(size=18,color="black"))+
    theme(plot.margin=unit(c(2,2,2,2), "lines"))
  #数据保存####
  names(prop)<- c("range","count")
  savexlsx(prop,paste0(savepath,proty,"_protein_sites.xlsx"),sheet = "sheet1")
  ggplotsave(plot = pro_plot,savepath=savepath,
             mapname = paste0(proty,"_protein_sites"),
             width = width,
             height = height,
             imagetype = imagetype,
             family=fontfamily)  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  #基本参数####
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 12, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 7.5, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  #此图参数####
  parser$add_argument("-if","--inputfile", default = NULL,help = "输入位点文件，默认为路径下的*Sites.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Statistics/", help = "输出结果路径，默认./Statistics/")
  args <- parser$parse_args()
  Prosite <- do.call(prosite,args = args)
}
