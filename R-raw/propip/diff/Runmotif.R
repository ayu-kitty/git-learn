#!/opt/conda/bin/Rscript

#' @export
premotif<-function(inputfile="差异表达矩阵.xlsx",savepath="./Motif/",group=group){
  createdir(savepath)
  diff <- readdata(inputfile,sheet = group)
  
  if("FlankingRegion" %in% colnames(diff)){
    sfas <- gsub(";.*","",diff$FlankingRegion) %>% na.omit()
  }else{
    sfas <- gsub(";.*","",diff$`Sequence window`)%>% na.omit()
  }
  
  sw <- nchar(sfas[1])
  write.table(sfas,paste0(savepath,"/",group,".txt"),sep="\t", 
              quote=FALSE, col.names=FALSE, row.names=FALSE)
  return(sw)
}

#' @export
motifmap<-function(compare=compare, imagetype = c("pdf", "png") , height = 4, width = 10, dpi = 300, fontfamily = "sans", ...){
  singlemotif<-function(i=i,momo=momo,number=number,compare=compare, imagetype = imagetype , height = height, width = width, dpi = dpi, fontfamily = fontfamily,...){
    pfm<-matrix(nrow=number,ncol=20)
    for(j in 1:number){
      mm<-as.vector(momo[(number+2)*i+j-number,])
      strsplit(mm,"\t")[[1]]->pfm[j,]
    }
    pfm<-apply(pfm,2,as.numeric)
    colnames(pfm)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    pfm<-t(pfm)
    colnames(pfm)<-as.character(1:number)
    name<-strsplit(as.vector(momo[(number+2)*i-(number+1),1])," ")[[1]][2]
    #sequence prob
    motprob<-ggseqlogo(pfm,method = "prob")+
      scale_y_continuous(expand = c(0,0))+
      theme(axis.ticks = element_line(color="black"))+
      theme(legend.text=element_text(size=12),legend.title = element_text(size=15),legend.position = "top")+
      theme(axis.text.x=element_text(size=12,color="black"),axis.text.y=element_text(size=12,color="black"),axis.title =element_text(size=15,color="black"))
    ggplotsave(plot = motprob,savepath=paste0(compare,"/","Probability"),
               mapname = name,
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily,
               ...)
    
    #sequence bits
    motbits<-ggseqlogo(pfm)+
      scale_y_continuous(expand = c(0,0))+
      theme(axis.ticks = element_line(color="black"))+
      theme(legend.text=element_text(size=12),legend.title = element_text(size=15),legend.position = "top")+
      theme(axis.text.x=element_text(size=12,color="black"),axis.text.y=element_text(size=12,color="black"),axis.title =element_text(size=15,color="black"))
    ggplotsave(plot = motbits,savepath=compare,
               mapname = name,
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily,
               ...)
    
    #sequence heatmap
  }
  if(!is.na(dir(compare,pattern="*.png")[1])){
    momo <- read.table(paste0(compare,"/momo.txt"),skip=6,sep="\n")
    number <- strsplit(momo[2,],"w= ")[[1]][2]%>%gsub(" .*","",.)%>%as.numeric()
    sapply(1:(nrow(momo)/(number+2)),function(i)singlemotif(i=i,momo=momo,number=number,compare=compare, imagetype = imagetype , height = height, width = width, dpi = dpi, fontfamily = fontfamily))
  }
}

#' run motif
#'
#' @param type 
#' @param number 
#' @param pval 
#' @param inputfile 
#' @param savepath 
#' @param group 
#' @param imagetype 
#' @param height 
#' @param width 
#' @param dpi 
#' @param fontfamily 
#' @param ... 
#' @export
runmotif <- function(number=NULL,pval=0.000001,
                     inputfile ="差异表达矩阵.xlsx",inputpath="./",savepath="./Motif/",group=NULL,
                     imagetype = c("pdf", "png") , height = 4, width = 10, dpi = 300, fontfamily = "sans", ...){
  pacman::p_load(openxlsx,ggseqlogo,ggplot2,dplyr)
  if(!file.exists(paste0(savepath,group,".txt"))){
    sw <- premotif(inputfile=paste0(inputpath,"/",inputfile),savepath=savepath,group=group)
  }else sw=31
  tmpdir <- paste0(savepath,"/",group)
  if(is.null(number)){
    dat=readdata(paste0(savepath,group,".txt"))
    number=ceiling(nrow(dat)/10)
  }
  if(number>=3){
    motifcommand <- "docker run --rm -v `pwd`:/mydata --user `id -u`:`id -g` -w /mydata abralab/memesuite:v5.1.1 momo motifx -oc "
    system(paste0(motifcommand,tmpdir," --verbosity 1 --width ",sw," --eliminate-repeats ",sw," --min-occurrences ",number," --score-threshold ",pval," ", paste0("'",savepath,group,".txt'")))
    ll <- motifmap(compare=tmpdir, imagetype = imagetype , height = height, width = width, dpi = dpi, fontfamily = fontfamily)
    print(paste0("~",group,"基序数量为：",length(dir(path = tmpdir,pattern = "*.png"))))
  }else{
    savetxt(data = paste0(group,"差异位点数量太少，不提供motif分析"),
            filename = paste0(savepath,"/说明.txt"),append = T)
    return()
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
  parser$add_argument("-he","--height",default = 4, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  parser$add_argument("-g","--group",default=NULL,help="比较组名称")
  parser$add_argument("-n","--number",default=NULL,help="筛选数量卡值，默认自动调整")
  parser$add_argument("-p","--pval",default=0.000001,help="筛选P值卡值，默认0.000001")
  parser$add_argument("-if","--inputfile",default="差异表达矩阵.xlsx",help="差异表格，默认为差异表达矩阵.xlsx")
  parser$add_argument("-ip","--inputpath",default = "./", help = "输入文件路径，默认./")
  parser$add_argument("-sp","--savepath",default = "./Motif/", help = "输出结果路径，默认./Motif/")
  
  args <- parser$parse_args()
  runmotif <- do.call(runmotif,args = args)
}
