#!/opt/conda/bin/Rscript

#' 空代堆积柱状图
#' 
#' @param areacolor 配色方案
#' @param filename 读入数据
#' @param imagetype 保存图片格式
#' @param savename 保存文件名
#' 
#' @export
stackedbarplot <- function(areacolor = SelectColors(palette = "customecol2"),
                       filename = "plotdata.xlsx",
                       imagetype = c("pdf","png"),
                       savename = "stackedbar_plot"){
  # 数据读入及长数据整理
  plotdata <- lmbio::readdata(filename=filename)
  plot <- melt(plotdata,id.vars="area")
  
  areaname <- colnames(plotdata)[-1]
  areacolor <- areacolor[1:length(areaname)]
  
  # 颜色整理
  ttest <- data.frame(areaname,areacolor)
  
  #按照柱子长度排序绘图
  orderclass<-ttest[order(ttest$areacolor),1] %>% as.factor()
  
  # 绘图
  pp=ggplot(plot, aes(x=area, y=value,fill=factor(variable,levels = orderclass)))+  
    geom_bar(stat = "identity",position="stack",width = 0.8,color="black")+theme_classic()+
    scale_fill_manual(values=setNames(areacolor,areaname))+
    theme(
      panel.grid=element_blank(),
      axis.text.x=element_text(size=15,angle =0,colour = "black"),
      plot.title=element_text(hjust=0.5,size=15),
      legend.title=element_blank(),
      panel.grid.major = element_blank(),
      legend.text = element_text(size=15),
      panel.grid.minor = element_blank(),
      axis.text.y=element_text(size=15,colour = "black"))+scale_y_continuous(expand = c(0,0))+
    labs(title ="",x="",y="")

  ggplotsave(pp,mapname = savename,imagetype = imagetype,height=8,width=12)
  
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-ac","--areacolor",default = SelectColors(palette = "customecol2"), help = "配色方案")
  parser$add_argument("-fn","--filename",default = "plotdata.xlsx", help = "读入数据")
  parser$add_argument("-i","--imagetype",default =  c("pdf","png"), help = "保存图片格式")
  parser$add_argument("-sn","--savename",default = "stackedbar_plot", help = "保存文件名")

  args <- parser$parse_args()
  
  writeinfo()
  
  mulargs <- do.call(what = stackedbarplot,args = args)
  
  writeinfo(endtime = T)
}

