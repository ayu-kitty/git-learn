#!/opt/conda/bin/Rscript

#' GO top30 plot
#' @param savepath 保存总富集路径
#' @param type 数据库类型缩写
#' @param incompare 比较组名称
#' @param filt 数据是否进行筛选
#' @param number 绘图top数量
#' @param imagetype 保存图片类型
#' @param height 高度
#' @param width 宽度
#' @param dpi 分辨率
#' @param fontfamily 字体类型
#' @param ... 
#' @export
enrichbarplot<-function(savepath = "./enrich/",indata="None", incompare = "A_B", type = "K",filt = "T",number =
                          10,imagetype = c("pdf","png") ,height = 8,width = 12,dpi=300,fontfamily="sans",...){
  pacman::p_load(dplyr,ggplot2,stringr,Hmisc)
  #colors = c("#FF7F0E","#1F77B4","#2CA03A","#FF7F0E","#1F77B4","#2CA03A","#D62728","#9068be","#8C564B","#D11250")
  colors <- c("#FF7F0E","#9068be","#2CA03A","#FF7F0E","#9068be","#2CA03A","#D62728","#17BECF","#D11250","#8C564B")
  gklevel1=c("biological_process","cellular_component","molecular_function",
             "Metabolism","Genetic Information Processing","Environmental Information Processing",
             "Cellular Processes","Organismal Systems","Human Diseases","Drug Development")
  gklevel1_color=data.frame("cateclass"=gklevel1,"Colors"=colors)
  if(indata=="None"){
    comenrichpath<-paste0(savepath,getenrichtype(type)$filn,"/",incompare,"/")
    enrichfiles<-paste0(comenrichpath,dir(path =comenrichpath, pattern = "enrichment-*"))
  }else{
    comenrichpath<-paste0(dirname(indata),"/")
    enrichfiles<-indata
  }
  for(enrichfile in enrichfiles){
    tudd<-ifelse(grepl("-Total",enrichfile),"Total",ifelse(grepl("-Up",enrichfile),"Up",ifelse(grepl("-Down",enrichfile),"Down","Total")))
    endat<-readdata(enrichfile)
    #### 列名报错提示
    correctcol <- c("Term","p-value","ListHits")
    findall <- grepl(paste0(colnames(endat),collapse = "|"),correctcol)
    if(any(!findall)){
      stop(paste0(paste0(correctcol[!findall],collapse = ","),"列名没有找到！请检查列名大小写"))
    }
    endat$`p-value`[endat$`p-value`==0]<-min(endat$`p-value`[endat$`p-value`!=0])
    if(length(grep("Category|Classification_level1",colnames(endat)))!=0){
      endat["cateclass"]<-endat[,grep("Category|Classification_level1",colnames(endat))]
    }
    if(filt == "T"){
      if(type == "G"){
        if(!("Category"%in%colnames(endat))){stop("Category列名没有找到！请检查列名大小写")}
        endata_filt<-rbind(filter(endat,ListHits!=1,cateclass=="molecular_function")[1:number,],
                      filter(endat,ListHits!=1,cateclass=="cellular_component")[1:number,],
                      filter(endat,ListHits!=1,cateclass=="biological_process")[1:number,])
        endata<-endata_filt %>% arrange(desc(cateclass), desc(`p-value`))
        endata_filt<-endata_filt[,-ncol(endata_filt)]
        endata<-left_join(endata,gklevel1_color,by="cateclass")
      }else if(type == "K"){
        if(!("Classification_level1"%in%colnames(endat))){stop("Classification_level1列名没有找到！请检查列名大小写")}
        endata_filt<-head(endat[order(endat$`p-value`),],number*2)
        endata<-endata_filt
        endata["catecount"]<-table(endata$cateclass)[endata$cateclass]
        endata<-endata %>% arrange(catecount,cateclass, desc(`p-value`))
        endata_filt<-endata_filt[,-ncol(endata_filt)]
        endata<-left_join(endata,gklevel1_color,by="cateclass")
      }else {
        endata_filt<-head(endat[order(endat$`p-value`),],number*2)
        endata<-endata_filt %>% arrange(desc(`p-value`))
      }
      na.omit(endata)->endata
      savexls(na.omit(endata_filt),paste0(comenrichpath,getenrichtype(type)$filn,".Bar.",tudd,".xls"))      
     }else {
      endata <- endat
      if(grepl(type,"G|K")){
        if(length(grep("Category|Classification_level1",colnames(endat)))==0){
          stop("Category或者Classification_level1列名没有找到！请检查列名大小写")
        }
        endata<-left_join(endata,gklevel1_color,by="cateclass")
      }
      endata<-as.data.frame(apply(endata,2,rev))
      na.omit(endata)->endata
    }
    if(nrow(endata)>2){
      endata$Term<-capitalize(as.character(endata$Term))
      endata$Term <- factor(endata$Term, levels=endata$Term)
      maxy<-max(-log10(as.numeric(endata$`p-value`)))
      if("cateclass" %in% colnames(endata)){
        p<-ggplot(data=endata, aes(y=Term, x=-log10(as.numeric(endata$`p-value`)), fill=cateclass)) +
          geom_bar(stat="identity",position=position_dodge(0.7),width=0.5, space=0.6,cex.main=3,cex.lab=2) +
          labs(y="", x="-log10(P-value)", title = ifelse(incompare!="",paste0(incompare,"(",tudd,")  \n Top ",getenrichtype(type)$filn," Terms"),
                                                   paste0("Top ",getenrichtype(type)$filn," Terms")),fill="") +
          xlim(0,1.08*maxy)+
          geom_text(aes(label=ListHits),hjust=-0.2, vjust=0.4,size=3) +
          theme_bw() + 
          theme(aspect.ratio=16/9)+
          theme(panel.grid=element_blank())+
          theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
          scale_fill_manual(values=rev(unique(endata$Colors)),
                            breaks =rev(unique(endata$cateclass)),
                            labels=rev(unique(endata$cateclass))) +
          theme(axis.text.y=element_text(hjust=1,size=12,color="black"),
                axis.text.x=element_text(size=12,color="black"),
                axis.title=element_text(size=14,color="black"),
                legend.text=element_text(size=14),legend.title=element_text(size=15))+
          scale_y_discrete(labels=function(x)ifelse(nchar(x)>60,paste0(substr(x,1,60),"..."),x))+
          theme(plot.title = element_text(hjust = 0.5, size=15))+
          theme(plot.margin=unit(c(2,2,2,2), "lines"))
      }else{
        p<-ggplot(data=endata, aes(y=Term, x=-log10(as.numeric(endata$`p-value`)))) +
          geom_bar(stat="identity",position=position_dodge(0.7),width=0.5, fill="#2CA03A",space=0.6,cex.main=3,cex.lab=2) +
          labs(y="", x="-log10 P-value", title = ifelse(incompare!="",paste0(incompare,"(",tudd,")  \n Top ",getenrichtype(type)$filn," Terms"),
                                                        paste0("Top ",getenrichtype(type)$filn," Terms"))) +
          xlim(0,1.08*maxy)+
          geom_text(aes(label=ListHits),hjust=-0.2, vjust=0.4,size=3) +
          theme_bw() + 
          theme(aspect.ratio=16/9)+
          theme(panel.grid=element_blank())+
          theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
          theme(axis.text.y=element_text(hjust=1,size=12,color="black"),
                axis.text.x=element_text(size=12,color="black"),
                axis.title=element_text(size=14,color="black"),
                legend.text=element_text(size=14),legend.title=element_text(size=15))+
          scale_y_discrete(labels=function(x)ifelse(nchar(x)>60,paste0(substr(x,1,60),"..."),x))+
          theme(plot.title = element_text(hjust = 0.5, size=15))+
          theme(plot.margin=unit(c(2,2,2,2), "lines"))
      }
      
      ggplotsave(plot = p,savepath=comenrichpath,
                 mapname = paste0(getenrichtype(type)$filn,".Bar.",tudd),
                 width = width,
                 height = height,
                 imagetype = imagetype,
                 family=fontfamily,
                 ...)
    }else{
      savetxt(data = "Term数量低于3，不提供富集柱状图",
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
  parser$add_argument("-wi","--width",default = 12, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 8, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-t","--type",type="character", default="K", help="enrich database class,such as G/K/W/R/I", metavar="character")
  parser$add_argument("-ic","--incompare",default = "A_B", help = "比较组名称，默认A_B")
  parser$add_argument("-id","--indata",default = "None", help = "绘图数据表格文件（包含路径），默认None不输入")
  parser$add_argument("-f","--filt",default = "T", type= "character",help = "数据是否需要筛选,默认T")
  parser$add_argument("-s","--savepath",default = "./enrich/", help = "富集结果保存路径,默认./enrich/")
  parser$add_argument("-n","--number",default = 10,help = "top number of term,默认10")
  args <- parser$parse_args()
  enbarplot <- do.call(enrichbarplot,args = args)
}

  