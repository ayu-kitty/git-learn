#!/opt/conda/bin/Rscript
#' 组内蛋白数鉴定分布图
#'
#' @param savepath 结果保存路径
#' @param inputpath 输入文件路径
#' @param inputfile 输入文件名称
#' @param sample sample表名称
#' @param imagetype 保存图片类型
#' @param height 图片高度
#' @param width 图片宽度
#' @param dpi 分辨率
#' @param fontfamily 字体样式
#' @param ... 
#'
#' @export
groupproplot <- function(savepath="./Identified/", inputpath="./",inputfile="Protein quantitation.xlsx",
                         sample = "Sample_Group.xlsx",imagetype=c("pdf","png"), 
                         classfile = "../classtype.xlsx",
                         height=9, width=10, dpi=300, fontfamily="sans", 
                         ...){
  pacman::p_load(ggplot2,dplyr,plotrix)
  createdir(filename = savepath)
  ##########数据整理##########
  data_pro <- readdata(paste0(inputpath,inputfile))
  if(length(unique(data_pro[,1]))!=nrow(data_pro)){
    protype="Sites"
  }else protype="Proteins"
  sample_infor <- readdata(paste0(inputpath,sample),sheet = 1)
  sample_infor <- sample_infor[which(sample_infor$Group!="删除"),]
  sg<-as.data.frame(table(sample_infor[,2]))
  if(length(as.character(sg[which(sg[,2]>1),1]))>0){
    sample_infor<-sample_infor[sample_infor[,2] %in% as.character(sg[which(sg[,2]>1),1]),]
    sample<-sample_infor[!(duplicated(sample_infor[,1])),]
    prod<-data_pro[,sample[,1]]
    rownames(sample)<-sample[,1]
    prod[prod == 0] <- NA
    gc<-as.data.frame(matrix(ncol=length(unique(sample_infor[,2])),nrow=nrow(prod)))
    gc<-sapply(1:length(unique(sample_infor[,2])), function(i){
      gs<-sample_infor[sample_infor[,2]==unique(sample_infor[,2])[i],1]
      sapply(1:nrow(prod), function(x){
        prod[x,gs]%>%as.numeric()%>%is.na()%>%sum()/length(gs)
      })
    }) %>% as.data.frame(gc)
    names(gc)<-unique(sample_infor[,2])
    gcount<-apply(gc, 2, function(x)length(which(x<=0.5)))
    data_group<-data.frame("group"=names(gcount),"number"=as.numeric(gcount))
    
    data_group[,1] <- factor(data_group[,1],levels = data_group[,1])
    sumax <- max(data_group[,2])*1.01
    ########绘图########
    group_plot <- ggplot(data_group,aes(group,number,fill=group))+
      geom_bar( stat="identity",width = 0.45)+
      labs(x="", y=paste0("Number of ",protype,"(Group)\n "),fill="") +
      ylim(0,sumax) +
      theme_bw() + 
      theme(panel.grid=element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
      geom_text(mapping = aes(label = number),vjust = -0.5,size=5)+
      scale_fill_manual(values=stylefun_group(classfile = classfile)[names(gcount)]) +
      theme(text=element_text(size=15),
            plot.margin = unit(rep(2,4),"lines"),
            legend.text=element_text(size=15),
            axis.text.x=element_text(angle=45,hjust = 1,size=15,color="black"),
            axis.text.y=element_text(size=15,color="black"),
            title = element_text(size = 18,color="black"))
    #########保存数据#########
    savexlsx(data_group,paste0(savepath,"Group_",protype,"_Number.xlsx"))
    ggplotsave(plot = group_plot,savepath=savepath,
               mapname = paste0("Group_",protype,"_Number"),
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily,
               ...)
  }
  
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
  parser$add_argument("-cf","--classfile",default = "../classtype.xlsx", help = "模板")
  
  # 此图参数
  parser$add_argument("-if","--inputfile", default = "Protein quantitation.xlsx",help = "定量结果文件名，默认为Protein quantitation.xlsx")
  parser$add_argument("-samp","--sample", default = "Sample_Group.xlsx",help = "包含样本信息的数据文件名，默认为Sample_Group.xlsx")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Identified/", help = "输出结果路径，默认./Identified/")
  args <- parser$parse_args()
  grouppro <- do.call(groupproplot,args = args)
}

