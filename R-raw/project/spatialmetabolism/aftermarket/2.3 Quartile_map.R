#!/opt/conda/bin/Rscript

#' 空代四分位图
#' 
#' @param areapath 文件路径
#' @param areacolor 配色方案
#' @param filename 读入数据
#' @param imagetype 保存图片格式
#' @param groupinfo 分组信息
#' @param needlist 代谢物名列
#' @param legend.position 图例位置
#' @param family 字体
#' @param imagetype 图片格式
#' @param comparisons 显著性标识比较组
#' @param classfile 分组文件
#' @param grouporder 分组排序
#' 
#' @export
Quartile_map <- function(filename="data.xlsx",
                         groupinfo = "groupdata.xlsx",
                         needlist = "Metabolites",
                         legend.position="right",
                         legend.text.size=8,
                         plot.title.size=14,
                         axis.text.size=10,
                         axis.title.size=13,
                         family="sans",
                         x.title="Group",
                         y.title="Expression",
                         height=6,
                         width=6,
                         dpi=300,
                         imagetype=c("jpg","pdf"),
                         grouporder=NULL,
                         comparisons=list(c("PC","ST"),c("PC","CT")),
                         classfile = "classtype.xlsx",
                         col = NULL){
  
  library(forcats)
  library(ggpubr)
  
  plotdata<-readdata(filename=filename,sheet=1)
  groupdata<-readdata(filename=groupinfo,sheet=1)
  
  if (is.null(grouporder)) {
    sample_group <- unique(groupdata$Group)
  } else {
    sample_group  <- unlist(strsplit(grouporder, split = ","))
  }
  
  if(is.null(comparisons)){
    group <- combn(sample_group, 2)
    comparisons <- list()
    for (i in 1:ncol(group)) comparisons[[i]] <- group[, i]
    
  }
  
  if(is.null(col)){
    col = stylefun_group(classfile = classfile,styletype = "fill")
  }
  
  for(i in 1:nrow(plotdata)){
    da1<-plotdata[i,]
    da2<-melt(da1,id.vars=needlist)
    da2<-merge(da2,groupdata,by.x="variable",by.y="sample")
    #计算最大/最小/均值/中位数
    da3<-da2 %>% group_by(Group) %>% mutate(upper=quantile(value,0.75),
                                            lower=quantile(value,0.25),
                                            mean=mean(value),
                                            median=median(value)
    )
    
    
    p<-ggplot(da3,aes(factor(Group,levels = sample_group),value,shape=Group))+
      geom_jitter(aes(color=Group),size=1,position=position_jitter(0.2))+
      scale_shape_manual(values=c(16,15,17))+
      scale_color_manual(values=col)+
      geom_errorbar(aes(ymin=lower,ymax=upper),width=0.1,size=0.5)+
      stat_summary(fun="mean",geom="crossbar",mapping=aes(ymin=..y..,ymax=..y..),
                   width=0.3,size=0.3)+labs(title=test$mz[i],x=x.title,y=y.title)+
      theme_classic()+
      theme(panel.grid=element_blank(),
            legend.position=legend.position,
            legend.text =element_text(size=legend.text.size,color="black",family=family),
            plot.title=element_text(size=plot.title.size,color="black",family=family,hjust=0.5),
            axis.text=element_text(size=axis.text.size,color="black",family=family),
            axis.title=element_text(size=axis.title.size,color="black",family=family),
            plot.margin = unit(c(1,1,1,1),"cm"))+
      
      stat_compare_means(method="t.test",hide.ns=F,
                         comparisons = comparisons,
                         label="p.signif",
                         bracker.size=0.8,
                         size=5
      )
    ggplotsave(p,mapname=test$mz[i],width=6,height=6,imagetype = imagetype,dpi=dpi)
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-ap","--areapath",default = SelectColors(palette = "customecol2"), help = "配色方案")
  parser$add_argument("-fn","--filename",default = "data.xlsx", help = "读入数据")
  parser$add_argument("-gi","--groupinfo",default = "groupdata.xlsx", help = "分组数据")
  parser$add_argument("-nl","--needlist",default = "Metabolites", help = "代谢物索引列")
  parser$add_argument("-i","--imagetype",default =  c("pdf","png"), help = "保存图片格式")
  parser$add_argument("-cp","--comparisons",default =  list(c("PC","ST"),c("PC","CT")), help = "显著性比较组")
  parser$add_argument("-cf","--classfile",default =  "classtype.xlsx", help = "分组数据")
  parser$add_argument("-fa","--family",default =  "sans", help = "字体")

  
  args <- parser$parse_args()
  
  writeinfo()
  
  mulargs <- do.call(what = Quartile_map,args = args)
  
  writeinfo(endtime = T)
}






