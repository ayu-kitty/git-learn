#!/opt/conda/bin/Rscript

#' @export
Pie_chart <- function(data = c("数据矩阵.xlsx","数据矩阵"),
                      mycolor2=c(
                        "#A6CEE3" ,"#1F78B4" ,"#B2DF8A" ,"#33A02C","#FB9A99" ,
                                  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"
                      ),
                      select = c("Super Class","Class","Sub Class"),
                      savepath = "./",
                      mapname = "piechart",
                      imagetype=c("jpg","pdf"),
                      height = 6,
                      width = 5,
                      ...){
  suppressMessages(library("Hmisc"))
  suppressMessages(library("dplyr"))
  
  datapie <- readdata(filename = data)
  datapie1 <- select(datapie,any_of(select))
  
  for(j in 1:dim(datapie1)[2]){
    datapie<-datapie1[,j][!datapie1[,j]=="Unclassified"] %>% as.data.frame()
    colnames(datapie)[1] <- colnames(datapie1)[j]
    datapie[,1] <- stringr::str_to_title(datapie[,1])
    datapie[,1]<-gsub(" And "," and ",datapie[,1])
    nclass=datapie %>%group_by(datapie[,1]) %>%count %>% ungroup
    colnames(nclass)[1]<-colnames(datapie[1])
    newdata =nclass[order(-nclass$n),]
    a<-data.frame(na.omit(newdata))
    if(length(a[,1])>9){
      a1<-a
      for(i in 1:length(a[,1])){
        a1[i,2]=a[i,2]
        if(a[i,2]< a[10,2]){
          a1[i,1]="Others"
        }else{a1[i,1]=a[i,1]}
      }
      
      a2<-data.frame()
      for(i in 1:9){
        a2[i,1]<-a1[i,1]
        a2[i,2]<-a1[i,2]
        a2[10,1]<-"Others"
        a2[10,2]<-sum(a1[10:length(a1[,1]),2])
      }
      colnames(a2)[1]<-colnames(datapie[1])
      colnames(a2)[2]<-"num"
      a2[,3]<-paste0(round(100*a2$num / sum(a2$num),2),"%")
      a3<-merge(nclass,a2,by=colnames(datapie[1]))
      a3<-a3[,c(1,3:4)]
      a3[10,]<-a2[10,]
      colnames(a3)[3]="Percent"
      colnames(a3)[2]="Count"
    }else{
      a1<-a
      a2<-data.frame()
      for(i in 1:length(a[,1])){
        a2[i,1]<-a1[i,1]
        a2[i,2]<-a1[i,2]
      }
      
      colnames(a2)[1]<-colnames(datapie[1])
      colnames(a2)[2]<-"num"
      a2[,3]<-paste0(round(100*a2$num / sum(a2$num),2),"%")
      a3<-merge(nclass,a2,by=colnames(datapie[1]))
      a3<-a3[,c(1,3:4)]
      colnames(a3)[3]="Percent"
      colnames(a3)[2]="Count"
    } 
    
    name1<-colnames(a3)[1]
    a3[,1]<- forcats::fct_inorder(a3[,1])
    a3$fraction<-a3$Count/sum(a3$Count)
    a3$ymax<-cumsum(a3$fraction)
    a3$ymin<-c(0,head(a3$ymax,n=-1))
    a3$label<-paste0(a3$Percent)
    a3$labelPosition<-(a3$ymax+a3$ymin)/2
    
    # filename1<-paste0(colnames(a3)[1],".jpg")
    # filename2<-paste0(colnames(a3)[1],".pdf")
    
    pp <- ggplot(a3,aes(ymax=ymax,ymin=ymin,xmax=2,xmin=1)) +
      geom_rect(aes(fill=a3[,1]),colour ="white")+
      scale_fill_manual(values = mycolor2)+
      coord_polar(theta = "y") +
      xlim(0,NA)+
      labs(x = "", y = "", title = "", fill = "") +
      theme_bw()+ 
      theme(axis.title = ggplot2::element_blank(),
            legend.text = element_text(size = 5,colour = "black"),
            legend.key.size = unit(0.8,"line"),
            axis.text =ggplot2::element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank())+
      annotate("text", x = 0, y = 0, label = name1, size = 5)+
      annotate("text", x = 2.15, y = a3$labelPosition, label =a3$Percent, size = 2)
    
    ggplotsave(plot = pp,
               savepath = savepath,
               mapname = paste0(mapname,"-",name1),
               imagetype = imagetype,
               width = width,
               height = height,
               ...)
    
    savexlsx1(data=a3[,1:3],
              filename = paste0(savepath,"/Pie_data.xlsx"),
              sheet = name1)
    
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_common_piechart <- map_autodraw$new(moudle = Pie_chart)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "饼图矩阵文件",required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf","html"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  result <- do.call(what = map_common_piechart ,args = args) 
}

#' 根据文件进行饼图可视化
#' 
#' @export
map_common_piechart <- map_autodraw$new(moudle = Pie_chart)$draw
