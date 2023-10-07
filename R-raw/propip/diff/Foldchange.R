#!/opt/conda/bin/Rscript

#' Title
#'
#' @param savepath 输出文件保存路径
#' @param upcolor 上调柱子颜色
#' @param downcolor 下调柱子颜色
#' @param inputfile 输入的差异表名称，默认'差异蛋白筛选结果.xlsx',可选参数'差异位点筛选结果.xlsx'
#' @param imagetype 保存图片类型
#' @param height 图片高度
#' @param width 图片宽度
#' @param dpi 分辨率
#' @param fontfamily 字体样式 
#' @param ... 
#'
#' @export
Foldchange<- function(savepath="./Foldchange/", inputpath="./", inputfile = '差异蛋白筛选结果.xlsx' ,
                      upcolor = "#f47d8c", downcolor = "#727fb5",
                      # upcolor = SelectColors(palette = "volcanocol",n = 5)[1],
                      # downcolor = SelectColors(palette = "volcanocol",n = 5)[5],
                      imagetype = c("pdf", "png") , height = 9, width = 10, dpi = 300, fontfamily = "sans", ...){
  pacman::p_load(openxlsx,readxl,stringr,ggplot2,reshape2)
  
  sheetn <- getsheetname(paste0(inputpath,"/",inputfile))
  com <- setdiff(sheetn,c(grep("原始",sheetn,value = T),grep("可信",sheetn,value = T)))
  #dataR为差异表格中上下调数量
  dataR <- data.frame()
  # d <- if (inputfile=="差异蛋白筛选结果.xlsx") {"Number Of Proteins"}else{
  #   if (inputfile=="差异位点筛选结果.xlsx") {"Number Of Sites"}else{"Number Of Peptides"}
  # }
  
  d <- "Number of Significant Differences"
  for (i in 1:length(com)) {
    if(com[i] == ""){
      for (j in 1:length(inputfile)) {
        diff <- readdata(paste0(inputpath,"/",inputfile[j]))
        comparename <- gsub(pattern = "/",replacement = "-vs-",x = diff[["args"]][["compare"]])
        diff <- get_diffana_data_2(data = diff)
        if("FoldChange"%in%colnames(diff)){
          dataR[comparename,1] <- comparename
          dataR[comparename,2] <- length(which(!as.numeric(diff$FoldChange)<1))
          dataR[comparename,3] <- length(which(as.numeric(diff$FoldChange)<1))
          dataR[comparename,4] <- nrow(diff)
        }else{
          dataR[comparename,1] <- comparename
          dataR[comparename,2] <- nrow(diff)
          dataR[comparename,3] <- NA
          dataR[comparename,4] <- nrow(diff)
        }
      }
    }else{
      diff <- readdata(paste0(inputpath,"/",inputfile),sheet=com[i])
      
      if("FoldChange"%in%colnames(diff)){
        dataR[sheetn[i],1] <- sheetn[i]
        dataR[sheetn[i],2] <- length(which(!as.numeric(diff$FoldChange)<1))
        dataR[sheetn[i],3] <- length(which(as.numeric(diff$FoldChange)<1))
        dataR[sheetn[i],4] <- nrow(diff)
      }else{
        dataR[sheetn[i],1] <- sheetn[i]
        dataR[sheetn[i],2] <- nrow(diff)
        dataR[sheetn[i],3] <- NA
        dataR[sheetn[i],4] <- nrow(diff)
      }
    }
  }
  colnames(dataR) <-  c("group","up","down","sum")
  # rownames(dataR) <- dataR$group
  dataR2 <- dataR
  colnames(dataR2) <-  c("Compare","Up","Down","Sum")
  savexlsx(data = dataR2,filename = paste0(savepath,"/差异统计.xlsx"))
  
  dataR <- dataR[,c("group","up","down")]
  com <- dataR$group
  if(length(com)>10){
    x <- ceiling(length(com)*0.1)
    y <- ceiling(length(com)/x)
	if(!x*y==(length(com))){
      com[(length(com)+1):(x*y)]<- NA
    }
    dd<- matrix(com, nrow = y, ncol = x)
  }else{
    dd <- matrix(com)
  }
  
  for(k in 1:ncol(dd)) {
    class <- dd[,k]
    class <- class[order(nchar(class))]
    data1 <- dataR[class,]
    df <-  melt(data1,id= "group")%>%na.omit()
    df$variable <- factor(df$variable,levels = c("up","down"))
    df$group <- factor(df$group,levels = c(class))
    if(length(com)==1) {
      p <- ggplot(df,aes(x=group,y=value,fill=variable))+
        geom_bar( stat="identity",position = position_dodge(0.6),width = 0.6)+
        labs(x = "", y = d)+
        geom_text(aes(label=value), size=5,
                  position = position_dodge(width = 0.6), vjust=-0.5)+
        scale_fill_manual(values = c(upcolor,downcolor))+
        theme(panel.grid = element_blank(),
              panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
              panel.background = element_blank(),
              axis.text.x = element_text(colour="black",hjust =1,angle = 45,size = ifelse(nchar(max(as.character(df$group)))>25,(ifelse(nchar(max(as.character(df$group)))>30,9,11)),15)),
              axis.text.y = element_text(colour="black",size = 16),
              axis.title.y = element_text(margin = margin(0,1,0,1,"cm"),size = 20),
              axis.ticks.y = element_line(color = "black",size = 1),
              legend.text = element_text(size = 14),
              plot.margin=unit(c(2,2,2,2), "lines"),
              aspect.ratio = 3/2)+
        guides(fill=guide_legend(title=""))
    }else{
      p <- ggplot(df,aes(x=group,y=value,fill=variable))+
        geom_bar( stat="identity",position = position_dodge(0.6),width = 0.6)+
        labs(x = "", y = d)+
        geom_text(aes(label=value), size=4,
                  position = position_dodge(width = 0.6), vjust=-0.5)+
        scale_fill_manual(values = c(upcolor,downcolor))+
        theme(panel.grid = element_blank(),
              panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
              panel.background = element_blank(),
              axis.text.x = element_text(colour="black",hjust =1,angle = 45,size = ifelse(nchar(max(as.character(df$group)))>25,(ifelse(nchar(max(as.character(df$group)))>30,9,11)),15)),
              axis.text.y = element_text(colour="black",size = 16),
              axis.title.y = element_text(margin = margin(0,1,0,1,"cm"),size = 20),
              axis.ticks.y = element_line(color = "black",size = 1),
              legend.text = element_text(size = 14),
              plot.margin=unit(c(2,2,2,2), "lines"))+
        guides(fill=guide_legend(title=""))
    }
    ggplotsave(plot = p,savepath=savepath,
               mapname = ifelse(ncol(dd)>1,paste0("foldchange_bars-",k),"foldchange_bars"),
               width = width,
               height = height,
               imagetype = imagetype,
               family = fontfamily)
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
  
  # 此图参数
  # parser$add_argument("-up","--upcolor",default = SelectColors(palette = "volcanocol",n = 5)[1], help = "上调柱子的颜色")
  # parser$add_argument("-down","--downcolor",default = SelectColors(palette = "volcanocol",n = 5)[5],help="下调柱子的颜色")
  parser$add_argument("-up","--upcolor",default = "#f47d8c", help = "上调柱子的颜色")
  parser$add_argument("-down","--downcolor",default = "#727fb5",help="下调柱子的颜色")
  parser$add_argument("-if","--inputfile",default="差异蛋白筛选结果.xlsx",help="差异表格，默认为差异蛋白筛选结果",nargs = "+")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "./Foldchange/", help = "输出结果路径，默认./Foldchange/")
  args <- parser$parse_args()
  Foldchange <- do.call(Foldchange,args = args)
}
