#!/opt/conda/bin/Rscript
#' @export
countlevel2<-function(pu){
  pu[!duplicated(pu[,2]),1:3]->level2
  level2[,"Count"]<-sapply(1:nrow(level2),function(x){
    pu[pu[,2]==level2[x,2],4]%>%strsplit(.,",")%>%unlist()%>%unique()%>%length})
  level2[,"Substances"]<-sapply(1:nrow(level2),function(x){
    pu[pu[,2]==level2[x,2],4]%>%strsplit(.,",")%>%unlist()%>%unique()%>%paste0(.,collapse =";")})
  return(level2)
}
#' KEGG level2上下调统计图
#' @param savepath 保存总富集路径
#' @param type 数据库类型缩写
#' @param incompare 比较组名称
#' @param filt 数据是否进行筛选
#' @param imagetype 保存图片类型
#' @param height 高度
#' @param width 宽度
#' @param dpi 分辨率
#' @param fontfamily 字体类型
#' @export
up_down<-function(savepath = "./enrich/", incompare = "A_B", filt = "T",  imagetype = c("pdf", "png") , height = 9, width = 14, dpi =
                    300, fontfamily = "sans", ...){
  pacman::p_load(dplyr,ggplot2,stringr,Hmisc)
  comenrichpath<-paste0(savepath,"KEGG/",incompare)
  upn<-dir(path =comenrichpath, pattern = "^enrichment.*-Up.xls")
  if(length(upn)==1){
    pathway_up<- read.table(paste0(comenrichpath,"/",upn),header=TRUE,sep="\t",check.names=F)
  }else pathway_up<-NULL
  down<-dir(path =comenrichpath, pattern = "^enrichment.*-Down.xls")
  if(length(down)==1){
    pathway_down<- read.table(paste0(comenrichpath,"/",down),header=TRUE,sep="\t",check.names=F)
  }else pathway_down<-NULL
  
  if(!is.null(pathway_up)){
    pu<-select(pathway_up,c("Classification_level1","Classification_level2","ListTotal","Substances"))
    lu<-cbind(countlevel2(pu),"Up")%>%as.data.frame()
    names(lu)[6]<-"Type"
  }else lu<-NULL
  if(!is.null(pathway_down)){
    pd<-select(pathway_down,c("Classification_level1","Classification_level2","ListTotal","Substances"))
    ld<-cbind(countlevel2(pd),"Down")%>%as.data.frame()
    names(ld)[6]<-"Type"
  }else ld<-NULL
  if(!is.null(lu) | !is.null(ld)){
    keggl<-rbind(lu,ld)%>%as.data.frame()
    if(nrow(keggl)>3){
      keggl[,"Percentage"]=keggl[,"Count"]/keggl[,"ListTotal"]*100
      arrange(keggl,Classification_level2)->kegglv
      write.table(kegglv,paste0(comenrichpath,"/Up_vs_Down.KEGG_Classification.xls"),sep="\t", row.names=F, quote=F)
      kegglv[,"mylabel"]=paste(kegglv$Classification_level1, kegglv$Classification_level2, sep="--")
      lbs<-which(table(kegglv[,"mylabel"])==1)%>% names()
      if(length(lbs)>0){
        kgna<-kegglv[kegglv[,"mylabel"] %in% lbs,]
        kgna[,"Count"]<-NA
        kgna[,"Percentage"]<-NA
        kgna[,"Type"]<-sapply(kgna[,"Type"], function(x)ifelse(x=="Up","Down","Up"))
        kegglv<-rbind(kegglv,kgna) %>% as.data.frame()
        arrange(kegglv,Classification_level2)->kegglv
      }
      
      kegglv$mylabel = factor(kegglv$mylabel,levels = unique(kegglv$mylabel),ordered = T)
      p<-ggplot(data=kegglv, aes(x=mylabel, y=Percentage, width=0.8, fill=Type, space=0)) +
        coord_flip() +
        geom_bar(stat="identity",position="dodge") +
        ylim(0,max(kegglv$Percentage)*1.1)+
        geom_text(aes(label=Count),position = position_dodge(width = 0.8), hjust=-0.5,size=4) + 
        scale_x_discrete(labels=function(x)ifelse(nchar(x)>60,paste0(substr(x,1,60),"..."),x))+
        guides(fill=guide_legend(reverse=TRUE))+
        scale_fill_manual(values=c('Down'='#727fb5','Up'='#f47d8c'))+
        labs(x="",y="Percent",fill="",title=paste0(incompare,"\n KEGG Pathway Classification"))+
        theme_bw()+
        theme(panel.grid=element_blank(),aspect.ratio=16/9)+
        theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
        theme(plot.title = element_text(hjust = 0.5, size=18)) +
        theme(plot.margin=unit(c(2,2,2,2), "lines"))+
        theme(axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black"),axis.title=element_text(size=18,color="black"),legend.text=element_text(size=14)) 
      
      ggplotsave(plot = p,savepath=comenrichpath,
                 mapname = paste0("/Up_vs_Down.KEGG_Classification"),
                 width = width,
                 height = height,
                 imagetype = imagetype,
                 family=fontfamily,
                 ...)
    }else{
      savetxt(data = "Term数量低于3，不提供level2富集上下调对比图",
              filename = paste0(comenrichpath,"/说明.txt"),append = T)
      return()
    }
  }
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("png","pdf"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 14, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 9, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-ic","--incompare",default = "A_B", help = "比较组")
  parser$add_argument("-f","--filt",default = "T", type= "character",help = "数据是否需要筛选,默认T")
  parser$add_argument("-s","--savepath",default = "./enrich/", help = "KEGG富集分析文件夹路径，默认./enrich/")
  args <- parser$parse_args()
  
  up_downplot <- do.call(up_down,args = args)
}

