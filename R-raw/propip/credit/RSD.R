#!/opt/conda/bin/Rscript
#' RSD分布图
#'
#' @param savepath 保存路径
#' @param inputpath 输入文件路径
#' @param samplefile Sample_Group.xlsx
#' @param inputfile 输入文件名称
#' @param imagetype 保存图片类型
#' @param height 高度
#' @param width 宽度
#' @param dpi 分辨率
#' @param fontfamily 字体类型 
#' @param classfile classtype.xlsx
#' @export
#'
prorsd <- function(savepath = "./Expression/RSD/",inputpath="./",# samplefile="Sample_Group.xlsx",
                   inputfile = "绘图数据.xlsx",imagetype = c("pdf","png") ,
                   height = 10,width = NULL,dpi=300,fontfamily="sans",
                   classfile = "../classtype.xlsx",
                   mapname = "RSD-boxplot",
                   ...){
  pacman::p_load(dplyr,ggplot2,scales,reshape2)
  
  #读取数据####
  # sample <- readdata(paste0(inputpath,"/",samplefile),sheet = 1)
  data <- readdata(paste0(inputpath,"/",inputfile))
  # rsdat <- select(data,all_of(sample[,1]))
  rsdat <- data[,-1,drop = F]
  rsdat <- 2^(rsdat)
  
  classdata <- readdata(paste0(dirname(classfile), "/classfile.yaml"))
  sample <- data.frame(sample = colnames(rsdat))
  for ( i in 1:dim(sample)[1]) {
    sample[i,"group"] <- names(classdata)[unlist(lapply(classdata, function(x,y){any(y %in% x)},y = sample[i,"sample"]))][1]
  }
  
  #数据处理####
  table <- as.data.frame(table(sample[,2]))#分组重复table
  
  if(max(table$Freq)==1){
    print("~本次实验所有组都无生物学重复，故不进行重复性检验！")
  }else{
    group<-table[table$Freq>1,1] %>% as.character()
    sample<-sample[sample[,2]%in%group,]
    
    rsddata <-sapply(1:nrow(rsdat), function(x){
      prod<-as.data.frame(cbind(t(rsdat[x,sample[,1]]),sample[,2]))
      sdd<-aggregate(as.numeric(prod[,1]),list(prod[,2]),sd)
      meand<-aggregate(as.numeric(prod[,1]),list(prod[,2]),mean)
      val<-sdd[2]$x/meand[2]$x
      names(val)<-sdd[[1]]
      return(val)
    }) %>% t() %>% as.data.frame()
    rsddata<-select(rsddata,unique(sample[,2]))
    rsddata_melt <- melt(rsddata,measure.vars = colnames(rsddata))
    colnames(rsddata_melt) <- c("Group","RSD")
    groupcol<-stylefun_group(classfile = classfile)[names(rsddata)]
    if(is.null(width)){
      width <- ifelse(ncol(rsddata)>80,6+2*(ceiling(ncol(rsddata)/20)),6+6*(ceiling(ncol(rsddata)/16)))
    }else{
      width <- width
    }
    #绘图####
	rsddata_melt$Group<-factor(rsddata_melt$Group,levels=unique(rsddata_melt$Group))
    p <- ggplot(rsddata_melt,aes(Group,RSD,color=Group))+
      stat_boxplot(geom = "errorbar",width=0.35)+
      geom_boxplot(aes(fill=Group),outlier.size=1)+
      geom_boxplot(alpha=0.5,outlier.size=1)+
      theme_bw() + 
      theme(panel.grid=element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1))+
      theme(plot.margin=unit(rep(2,4),'lines'))+
      theme(axis.text.x = element_text(size = 15,angle = 45,hjust = 1,color="black"),
            axis.text.y = element_text(size = 15,color="black"),
            axis.title = element_text(size = 18,color="black"),
            legend.text=element_text(size=15))+
      scale_color_manual(values=groupcol)+
      scale_fill_manual(values=groupcol)+
      guides(color=guide_legend(title=NULL),fill="none")+
      labs(x="", y="RSD",title="",fill="")
    #保存数据####
    if(colnames(data)[1] == "Accession"){
      result_data <- cbind(data["Accession"],rsddata)
      savexlsx(result_data,paste0(savepath,"/RSD-boxplot.xlsx"))
    }
    ggplotsave(plot = p,savepath=savepath,
               mapname = mapname,
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily)
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("png","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf"))
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-cf","--classfile",default = "../classtype.xlsx", help = "模板")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-if","--inputfile",type="character", default="绘图数据.xlsx", help="需要分析的表格，默认绘图数据.xlsx", metavar="character")
  parser$add_argument("-ip","--inputpath", default="./", help="分析文件路径，默认当前路径")
  # parser$add_argument("-samp","--samplefile",type="character", default="Sample_Group.xlsx", help="分组信息表格，默认Sample_Group.xlsx", metavar="character")
  parser$add_argument("-sp","--savepath",default = "./Expression/RSD/", help = "结果保存路径,默认./Expression/RSD/")
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  Rsd <- do.call(prorsd,args = args)
}
