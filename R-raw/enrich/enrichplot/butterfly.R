#!/opt/conda/bin/Rscript

#' 蝴蝶图
#'
#' @param savepath 保存路径
#' @param type 数据库类型
#' @param incompare 比较组名称
#' @param filt 是否进行表格筛选
#' @param imagetype 保存图片类型
#' @param height 高度
#' @param width 宽度
#' @param dpi 分辨率
#' @param fontfamily 字体类型 
#' @param ... 
#' @export
butterfly<-function(savepath = "./enrich/", type = "K", incompare = "A_B", filt = "T", imagetype = c("pdf", "png") , height = 8, width = 14, dpi =
                      300, fontfamily = "sans", ...){
  pacman::p_load(dplyr,ggplot2,stringr,Hmisc)
  
  path<-paste0(savepath,getenrichtype(type)$filn,"/",incompare)
  upf<-dir(path =path, pattern = "^enrichment.*-Up.xls")
  downf<-dir(path =path, pattern = "^enrichment.*-Down.xls")
  if(length(upf)==1 & length(downf)==1){
    data_up<- read.delim(paste0(path,"/",upf), header=T, sep="\t", quote="",check.names=FALSE)
    data_down<- read.delim(paste0(path,"/",downf), header=T, sep="\t", quote="",check.names=FALSE)
    if(filt=="T"){
      data_up<- data_up[which(data_up$ListHits>1),]
      data_down<- data_down[which(data_down$ListHits>1),]
    }
    if(nrow(data_down)>2 & nrow(data_up)>2){
      down<- head(data_down[order(data_down$`p-value`),],10)[,c("Term","ListHits","ListTotal","p-value")]##选取分别以down和up的pvalue排序前10的条目做数据框绘图
      down[,"Percent"]<-down[,"ListHits"]/down[,"ListTotal"]
      down[,"Type"]<- "Down"
      up<- head(data_up[order(data_up$`p-value`), ], 10)[,c("Term","ListHits","ListTotal","p-value")]
      up[,"Percent"]<-up[,"ListHits"]/up[,"ListTotal"]
      up[,"Type"]<- "Up"
      data<- rbind(down,up)%>%na.omit()
      write.table(data,paste0(path,"/Up_Down.Comparison.xls"),sep="\t", row.names=F, quote=F)
      data <- data %>%
        mutate(List = ifelse(Type == "Down", -Percent, Percent)) %>%arrange(Type,List)
      data$Term<-capitalize(as.character(data$Term))
      data[,"Labels"]<-sapply(as.character(data$Term),function(x){ifelse(nchar(x)>50,paste0(substr(x,1,50),"..."),x)})
      data$Term = factor(data$Term,levels = unique(data$Term),ordered = T)
      tmp = with(data, labeling::extended(range(-Percent)[1]-0.1, range(Percent)[2]+0.1, m = 5,only.loose = TRUE))
      lm = tmp[c(1,length(tmp))]
      p<- ggplot(data, aes(y=Term, x=List,fill=Type))+
        geom_bar(stat='identity',width=.7)+
        scale_fill_manual(values = c('Up'='#f47d8c','Down'='#727fb5'))+
        #scale_fill_manual(values = c('Up'='#DE897B','Down'='#9BCCE3'))+
        scale_x_continuous(limits = lm,breaks = tmp,labels = abs(tmp))+
        labs(y="", x="Percent", title = paste0(incompare,"\n ",getenrichtype(type)$filn," Up-Down Comparison Terms"))+
        geom_text(data = subset(data, List < 0),aes(y=Term,x=(tmp[length(tmp)]-tmp[1])/100, label=Labels),size = 5,hjust = 0)+
        geom_text(data = subset(data[!duplicated(data$Term),],List > 0),aes(y=Term,x=-(tmp[length(tmp)]-tmp[1])/100,label= Labels),size = 5, hjust = 1)+
        theme(aspect.ratio=9/16,
              panel.grid =element_blank(),#去除网格线
              panel.background = element_blank(),#背景颜色
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(), #去除y轴
              axis.text.y = element_blank(),
              axis.line.x = element_line(color = "black"),
              legend.title = element_blank(),
              legend.text = element_text(size = 14),
              axis.text = element_text(size=14,color="black"),
              axis.title = element_text(size=15,color="black"),
              plot.title = element_text(size=16,hjust=0.5,color="black"),
              plot.margin =unit(c(2,2,2,2), "lines"))
      
      ggplotsave(plot = p,savepath=path,
                 mapname = paste0("/Up_Down.Comparison"),
                 width = width,
                 height = height,
                 imagetype = imagetype,
                 family=fontfamily,
                 ...)
    }else{
      savetxt(data = "Term数量低于2，不提供富集上下调对比图",
              filename = paste0(path,"/说明.txt"),append = T)
      return()
    }
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
  parser$add_argument("-t","--type",type="character", default="K", help="enrich database class,such as G/K/W/R/I", metavar="character")
  parser$add_argument("-ic","--incompare",default = "A_B", help = "比较组")
  parser$add_argument("-f","--filt",default = "T", help = "是否需要筛选")
  parser$add_argument("-s","--savepath",default = "./enrich/", help = "KEGG富集分析文件夹路径，默认./enrich/")
  args <- parser$parse_args()
  butterflyplot <- do.call(butterfly,args = args)
}

