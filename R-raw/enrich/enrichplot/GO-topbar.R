#!/opt/conda/bin/Rscript

#' top20 GO barplot
#'
#' @param savepath 保存总富集路径
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
gotopbar<-function(savepath = "./enrich/", incompare = "A_B", filt = "T",number =
                     20,imagetype = c("pdf","png") ,height = 8,width = 12,dpi=300,fontfamily="sans",...){
  pacman::p_load(dplyr,ggplot2,stringr,Hmisc)
  
  BCM <- c("biological_process","cellular_component","molecular_function")
  fillcol <- c("#FF7F0E","#9068be","#2CA03A")
  comenrichpath<-paste0(savepath,"GO/",incompare)
  enrichfiles<-dir(path =comenrichpath, pattern = "enrichment-*")
  createdir(paste0(comenrichpath,"/GO-top",number))
  for(i in 1:length(enrichfiles)){
    ll<-strsplit(enrichfiles[i],"-")[[1]]
	ll[length(ll)]%>%strsplit(.,".xls")%>%unlist->tudd
    got<-read.delim(paste0(comenrichpath,"/",enrichfiles[i]), sep="\t", header=T, quote="",check.names=FALSE)
    for(k in 1:3){
      if(filt=="T"){
        go<-filter(got,ListHits!=1,Category==BCM[k])[1:number,]
        na.omit(go)->go		
        go<-go %>% arrange(Category, desc(ListHits))
        write.table(go,paste0(comenrichpath,"/GO-top",number,"/",BCM[k],"-",tudd,".xls"),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
      }else go<-got
      
      na.omit(go)->go
      go<-go %>% arrange(Category, desc(ListHits))
      if(nrow(go)>3){
        go$Term<-capitalize(as.character(go$Term))
        plotdata_top <- go[order(go$ListHits),]
        plotdata_top$Term <- factor(plotdata_top$Term, levels=plotdata_top$Term)
        p <- ggplot(plotdata_top,aes(x=ListHits,y=Term))+
          geom_bar(stat="identity",position=position_dodge(0.7),width=0.5, fill=fillcol[k],col=fillcol[k],space=0.6,cex.main=3,cex.lab=2)+
          theme_bw()+
          theme(panel.grid=element_blank())+
          theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
          labs(x="Protein number", y="",title=paste0(incompare," (",tudd,")","\n Top ",number," GO Terms (",BCM[k],")"))+
          theme(axis.text = element_text(size=15,color="black"),
                axis.title = element_text(size=18,color="black"),plot.title = element_text(size=18,hjust=0.12,color="black")) +
          theme(legend.position = "none")+
          theme(aspect.ratio=4/3)+
          theme(plot.margin=unit(c(2,2,2,2), "lines"))+
          scale_y_discrete(labels=function(x)ifelse(nchar(x)>60,paste0(substr(x,1,60),"..."),x))
        
        ggplotsave(plot = p,savepath=comenrichpath,
                   mapname = paste0("GO-top",number,"/",BCM[k],"-",tudd),
                   width = width,
                   height = height,
                   imagetype = imagetype,
                   family=fontfamily,
                   ...)
        
      }else{
        savetxt(data = "Term数量低于3，不提供GO category绘图",
                filename = paste0(comenrichpath,"/GO-top",number,"/","说明.txt"),append = T)
        return()
      }
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
  
  parser$add_argument("-ic","--incompare",default = "A_B", help = "比较组名称，默认A_B")
  parser$add_argument("-f","--filt",default = "T", type= "character",help = "数据是否需要筛选,默认T")
  parser$add_argument("-s","--savepath",default = "./enrich/", help = "富集结果保存路径,默认./enrich/")
  parser$add_argument("-n","--number",default = 20,help = "top number of term,默认20")
  args <- parser$parse_args()
  
  gotopbarplot <- do.call(gotopbar,args = args)
}