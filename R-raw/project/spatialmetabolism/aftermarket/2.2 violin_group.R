#!/opt/conda/bin/Rscript

#' 空代分组小提琴图
#' 
#' @param areapath 文件路径
#' @param areacolor 配色方案
#' @param filename 读入数据
#' @param imagetype 保存图片格式
#' @param groupinfo 分组信息
#' 
#' @export
violin_group <- function(areapath = "./",
                         filelist = list.files(path=areapath,pattern = ".txt"),
                         groupinfo = "groupdata.xlsx",
                         areacolor = c("#00468BFF","#ED0000FF","#ADB6B6FF"),
                         imagetype = c("jpg", "pdf")){
  
  groupsheet<-gsub(".txt","",filelist)
  groupinfo<-"groupdata.xlsx"
  groupdata <- readdata(groupinfo)
  areaname <- unique(groupdata$Group)
  for(i in 1:length(filelist)){
    #读取选区数据和分组方式，合并数据
    areadata<-lmbio::readdata(filename=paste0(areapath,filelist[i]))
    groupdata<-lmbio::readdata(filename=groupinfo,sheet=groupsheet[i])
    newdata<-melt(areadata,id.vars="mz")
    plotdata<-merge(newdata,groupdata,by.x="variable",by.y="sample")
    
    #按照每一个mz进行绘图
    for(mz in unique(plotdata$mz)){
      plotdata1<-plotdata[plotdata$mz==mz,]
      pp=ggviolin(plotdata1,x="area",y="value",fill="group",legend = "right",legend.title = "Group")+
        labs(title="",y="Expression",x="Group")+
        scale_y_continuous(expand = expansion(mult = .1))+
        scale_fill_manual(values=setNames(areacolor,areaname))+
        theme(
          plot.title=element_text(hjust=0.5,size=15),
          legend.text =element_text(size=9,color="black",angle=0),
          axis.text.x=element_text(size=11,color="black",angle=0),
          axis.text.y=element_text(size=11,color="black",angle=0),
          axis.title.y=element_text(size=13,color="black",angle=90,vjust=1.5),
          axis.title.x=element_text(size=13,color="black",angle=0,vjust=0.5),
          plot.margin = unit(c(1,1,1,1), "cm"))
      ggplotsave(pp,mapname=paste0("violin-",mz),height=8,width=12,imagetype = imagetype)
    }
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-ap","--areapath",default = "./", help = "文件路径")
  parser$add_argument("-fl","--filelist",default = "plotdata.xlsx", help = "读入数据")
  parser$add_argument("-gi","--groupinfo",default =  c("pdf","png"), help = "保存图片格式")
  parser$add_argument("-ac","--areacolor",default = "stackedbar_plot", help = "保存文件名")
  parser$add_argument("-i","--imagetype",default =  c("pdf","png"), help = "保存图片格式")
  
  args <- parser$parse_args()
  
  writeinfo()
  
  mulargs <- do.call(what = violin_group,args = args)
  
  writeinfo(endtime = T)
}


