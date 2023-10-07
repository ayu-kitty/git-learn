#!/opt/conda/bin/Rscript
#' 样本相关性图
#' 
#' @param savepath 保存路径
#' @param inputpath 输入文件路径
#' @param inputfile 输入文件名称
#' @param classfile classtype.xlsx
#' @param imagetype 保存图片类型
#' @param height 高度
#' @param width 宽度
#' @param dpi 分辨率
#' @param fontfamily 字体类型 
#'
#' @export
#'
samplecorr <- function(savepath = "./Expression/Samplecorrplot/",inputpath="./",
                       inputfile="绘图数据.xlsx",imagetype = c("png", "pdf"),
                       height = 18,width = 20,dpi=300,fontfamily="sans",
                       classfile = "../classtype.xlsx",
                       meta = F,
                       mapname = "Samplecorrplot",
                       ...){
  pacman::p_load(dplyr,corrplot,ggplot2,reshape,openxlsx)
  cols=c("#00468BFF","#0099B4FF","cornsilk3","#f19b78","#E64B35FF")
  #读取分组信息
  data <- readdata(paste0(inputpath,inputfile))[,-1]
  if(ncol(data)>2){
    if(meta){
      classdata <- readdata(paste0(dirname(classfile), "/classfile.yaml"))
      data <- data[,classdata[["QC"]],drop = F]
      samplecol <- rep("black",dim(data)[2])
    }else{
      samplecol <- stylefun_group(sample = T,classfile = classfile)[names(data)]
    }
    
    createdir(filename = savepath)
    ##参数设置
    sample_number <- c(2:200)
    number0 <- seq(48,0.5,length.out = 199)
    number0[6:9] <- number0[6:9]-8
    number0[4:5] <- number0[4:5]-15
    number0[2:3] <- number0[2:3]-20
    number0[1] <- number0[1]-45
    tl.cex <- number0[which(length(samplecol)==sample_number)]/length(samplecol)
    clcex=2
    
    #空值最小值填充
    if(is.na(sum(data))){
      data[is.na(data)]<-apply(data,2,function(x){min(na.omit(x))})%>%min()
    }
    corr_data <- as.matrix(cor(data,method = "pearson"))
    res1 <- cor.mtest(corr_data,conf.level = .95)
    if (length(sample)>100 & length(sample)<=200) {
      plotfile(savepath = savepath,mapname = mapname,imagetype = imagetype, 
               width = width, height = height,family=fontfamily)
      corrplot::corrplot(corr_data,type="upper",tl.pos="lt",order="original",
                         insig = "label_sig",
                         sig.level = c(.001,.01,.05),
                         col=cols,addgrid.col = 'grey33',
                         method = "circle",
                         tl.cex = tl.cex/(0.2),cl.cex = clcex,
                         tl.col = samplecol,mar = c(6, 6, 6, 6))
      corrplot::corrplot(corr_data,add=TRUE,type="lower", method="number",
                         order="original",diag=F,
                         col=cols,addgrid.col = 'grey33',
                         number.cex = tl.cex/(0.5),cl.cex = clcex,
                         cl.pos = "n",tl.pos = "n",mar = c(6, 6, 6, 6))
      plotsave()
    }else if (length(sample)>200) {
      plotfile(savepath = savepath,mapname = mapname,imagetype = imagetype, 
               width = width, height = height,family=fontfamily)
      corrplot::corrplot(corr_data,type="upper",tl.pos="lt",order="original",
                         insig = "label_sig",
                         sig.level = c(.001,.01,.05),
                         col=cols,addgrid.col = 'grey33',
                         method = "circle",
                         tl.cex = 0.01,cl.cex = clcex,
                         tl.col = samplecol,mar = c(6, 6, 6, 6))
      corrplot::corrplot(corr_data,add=TRUE,type="lower", method="number",
                         order="original",diag=F,
                         col=cols,addgrid.col = 'grey33',
                         number.cex = 0.01,cl.cex = clcex,
                         cl.pos = "n",tl.pos = "n",mar = c(6, 6, 6, 6))
      plotsave()
    }else{
      plotfile(savepath = savepath,mapname = mapname,imagetype = imagetype, 
               width = width, height = height,family=fontfamily)
      corrplot::corrplot(corr_data,type="upper",tl.pos="lt",order="original",
                         p.mat = res1$p,
                         insig = "label_sig",
                         sig.level = c(.001,.01,.05),
                         pch.col = "white",
                         pch.cex = tl.cex/(6/5),
                         col=cols,addgrid.col = 'grey33',
                         method = "circle",
                         tl.cex = tl.cex/(5/4),cl.cex =clcex,
                         tl.col = samplecol,mar = c(6, 6, 6, 6))
      corrplot::corrplot(corr_data,add=TRUE,type="lower", method="number",
                         order="original",diag=F,
                         col=cols,addgrid.col = 'grey33',
                         number.cex = tl.cex/(6/4),cl.cex = clcex,
                         cl.pos = "n",tl.pos = "n",mar = c(6, 6, 6, 6))
      plotsave()
    }
    wb <- createWorkbook()
    wb <- addsheet3(data = corr_data,wb = wb,sheet = "相关性系数")
    wb <- addsheet3(data = res1$p,wb = wb,sheet = "相关性检验")
    savewb(wb = wb, filename = paste0(savepath,"/",mapname,".xlsx")) 
  }
  
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 20, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 18, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf"))
  parser$add_argument("-cf","--classfile",default = "../classtype.xlsx", help = "模板")
  
  # 此图参数
  parser$add_argument("-mn","--mapname", default = "Samplecorrplot", help = "保存文件名")
  parser$add_argument("-if","--inputfile",type="character", default="绘图数据.xlsx", help="需要分析的表格，默认绘图数据.xlsx", metavar="character")
  parser$add_argument("-ip","--inputpath", default="./", help="分析文件路径，默认当前路径")
  parser$add_argument("-s","--savepath",default = "./Expression/Samplecorrplot/", help = "结果保存路径,默认./Expression/Samplecorrplot/")
  parser$add_argument("-m","--meta", default = F, action = "store_true", help = "是否是代谢数据")
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  Samplecorr <- do.call(samplecorr,args = args)
  
}
