#!/opt/conda/bin/Rscript
#' 联合分析火山图
#'
#' @param savepath 结果保存路径
#' @param inputpath 输入数据路径
#' @param backpath 注释数据路径
#' @param vesion 输入数据名称
#' @param imagetype 图片类型
#' @param dpi 图片分辨率
#' @param height 图片高度
#' @param width 图片宽度
#' @param proFC 蛋白组FC
#' @param phosFC 修饰组FC
#' @param pvalue P-value值
#' @param fontfamily 字体样式
#' @param ... 
#' @export
unionvolplot <- function(savepath="Scatterplot/", inputpath="rawdata/",backpath="background/",vesion="差异表达矩阵(未筛选)",
                         imagetype=c("pdf","png"), dpi=300, height=8,width=10,
                         proFC = NULL,phosFC= NULL,pvalue = NULL,
                         fontfamily="sans", ...){
  pacman::p_load(openxlsx,readxl,stringr,ggplot2,dplyr)
  if(file.exists(paste0(inputpath,"ratio.xlsx"))){
    ratio <- readdata(paste0(inputpath,"ratio.xlsx"))
    proFC <- ratio[1,"FC_pro"]
    phosFC <- ratio[1,"FC_ptm"]
    pvalue <- ratio[1,"P-value_ptm"]
  }else print("~！！！！rawdata路径下不存在ratio.xlsx文件，请输入FC和p值参数！！！")
  cg<-excel_sheets(paste0(inputpath,vesion,"-pro.xlsx"))
  anno <- readdata(paste0(backpath,"Annotation_Data.xlsx"))[,-1]
  samp <- readdata(paste0(inputpath,"Sample_Group.xlsx"),sheet=1)
  names(anno)[1] <- "Accession"
  wb <- openxlsx::createWorkbook()
  for (i in 1:length(cg)) {
    #读取数据
    prodata <- readdata(paste0(inputpath,vesion,"-pro.xlsx"),sheet = cg[i])
    names(prodata) <- gsub("p-value","P-value",names(prodata))
    names(prodata) <- gsub("FoldChange","FC",names(prodata))
    proteindata <- select(prodata,one_of(c("Accession","Gene Name","P-value","FC")))
    proteindata$Regulation <- "None"
    
    proteindata$Regulation[proteindata$FC > 1] <- "Up"
    proteindata$Regulation[proteindata$FC < 1] <- "Down"
    names(proteindata)[-1:-2] <- paste(names(proteindata)[-1:-2],"protein",sep = "_")
    site <- readdata(paste0(inputpath,vesion,"-ptm.xlsx"),sheet = cg[i])

    ptmn<-ptma(colnames(site))
    names(site) <- gsub("p-value","P-value",names(site))
    names(site) <- gsub("FoldChange","FC",names(site))
    phosdata <- select(site,one_of(c("Accession","Positions within proteins","Amino acid","P-value","FC")))
    phosdata[,"prosite"]<-gsub(";.*","",phosdata$`Positions within proteins`)%>%paste0(phosdata$Accession,":",phosdata$`Amino acid`,.)
    phosda<-select(phosdata,one_of(c("Accession","prosite","P-value","FC")))
    phosda$Regulation <- "None"
    if("P-value" %in% colnames(phosda)){
      phosda$Regulation[phosda$FC >= phosFC & phosda$`P-value` < pvalue] <- "Up"
      phosda$Regulation[phosda$FC <= 1/phosFC & phosda$`P-value` < pvalue] <- "Down"
    }else{
      phosda$Regulation[phosda$FC >= phosFC] <- "Up"
      phosda$Regulation[phosda$FC <= 1/phosFC] <- "Down"
    }
    names(phosda)[-1:-2] <- paste(names(phosda)[-1:-2],paste0(ptmn,"protein"),sep = "_")
    sitefc <- paste0("FC_",ptmn,"protein")
    alldata <- merge(x = proteindata,y = phosda)
    alldata$FC_protein <- as.numeric(alldata$FC_protein)
    alldata[,sitefc] <- as.numeric(alldata[,sitefc])
    alldata[alldata$FC_protein == 0,"FC_protein"] <- 0.0078125
    alldata[alldata$FC_protein == Inf,"FC_protein"] <- 128
    alldata[alldata[,sitefc] == 0,sitefc] <- 0.0078125
    alldata[alldata[,sitefc] == Inf,sitefc] <- 128
    if("P-value_protein" %in% colnames(alldata)){
      alldata <- alldata[,c(6,1,2,8,7,9,4,3,5)]
      alldata[,"class"] <- "not-significant"
      alldata$class[alldata[,6]=="Up" & (alldata[,9]=="Down")] <- "change-up"
      alldata[alldata[,6]=="Up" & alldata[,9]=="Up" & alldata[,4]-alldata[,7] >0 ,"class"] <- "more-up"
      
      alldata$class[alldata[,6]=="Down" & (alldata[,9]=="Up")] <- "change-down"
      alldata[alldata[,6]=="Down" & alldata[,9]=="Down" & alldata[,4]-alldata[,7] <0 ,"class"] <- "more-down"
    }else{
      alldata <- alldata[,c(5,1,2,6,7,3,4)]
      alldata[,"class"] <- "not-significant"
      alldata$class[alldata[,5]=="Up" & (alldata[,7]=="Down")] <- "change-up"
      alldata[alldata[,5]=="Up" & alldata[,7]=="Up" & alldata[,4]-alldata[,6] >0 ,"class"] <- "more-up"
      
      alldata$class[alldata[,5]=="Down" & (alldata[,7]=="Up")] <- "change-down"
      alldata[alldata[,5]=="Down" & alldata[,7]=="Down" & alldata[,4]-alldata[,6] <0 ,"class"] <- "more-down"
    }
    vol <- alldata
    vol$FC_protein <- log2(vol$FC_protein)
    vol[,sitefc] <- log2(vol[,sitefc])
    vol$class <- factor(x = vol$class,levels = c("not-significant","change-up","more-up","change-down","more-down"))
    vol <- arrange(vol,class)
    
    x <- min(c(max(abs(c(quantile(vol[,"FC_protein"], probs = c(0.9999))*1.1,quantile(vol[,"FC_protein"], probs = c(0.0001))*1.1))),15))
    y <- min(c(max(abs(c(quantile(vol[,sitefc], probs = c(0.9999))*1.1,quantile(vol[,sitefc], probs = c(0.0001))*1.1))),15))
    a <- ceiling(max(x,y))
    ##########绘图##########
    ldata <- data.frame(mx=seq(-a,a,by=0.1),my=seq(-a,a,by=0.1))
    pp <- ggplot(data = vol,mapping = aes(x = FC_protein,y = vol[,sitefc],color=class))+
      geom_point()+
      geom_line(data = ldata,mapping = aes(x=mx,y=my),
                color = "black", size = 0.5) +
      #geom_abline(intercept=0,slope=1,color = "black", size = 0.5) +
      scale_x_continuous(limits = c(-x,x))+
      scale_y_continuous(limits = c(-y,y))+
      scale_color_manual(values=c("change-up" = "#FEBFCA","more-up" = "#FF4D40","not-significant" = "#D2DAE2","change-down" = "#86CDF9","more-down" = "#4682B4"),
                         name="Regulation Type",
                         breaks = c("change-up","more-up","change-down","more-down","not-significant"),
                         labels = c(paste0("protein-down&",ptmn,"-up"),
                                    paste0(ptmn,"-more up"),
                                    paste0("protein-up&",ptmn,"-down"),
                                    paste0(ptmn,"-more down"),
                                    paste0("not-significant")))+
      guides(size = guide_legend(order = 1))+
      labs(x="log2(Protein FC)",y=paste0("log2(",ptmn,"protein FC)"),title = paste0("Volcano Plot (",cg[i],")"))+
      geom_vline(xintercept = c(log2(proFC)),linetype="dashed",size=0.5,colour="black")+
      geom_vline(xintercept = c(-log2(proFC)),linetype="dashed",size=0.5,colour="black")+
      geom_hline(yintercept = c(log2(phosFC)),linetype="dashed",size=0.5,colour="black")+
      geom_hline(yintercept = c(-log2(phosFC)),linetype="dashed",size=0.5,colour="black")+
      theme(panel.grid=element_blank(),panel.background = element_rect(fill = NA),
            panel.border = element_rect(fill=NA,color="black",size = 1,linetype="solid"),
            axis.text.x=element_text(size=15,color="black"),
            axis.text.y=element_text(size=15,color="black"),
            plot.margin = ggplot2::unit(c(0.3,0.3,0.3,0.3),"in"),
            plot.title = element_text(hjust = 0.5,size = 18),
            legend.background = element_rect(),
            legend.key = element_blank(),
            text=element_text(size=15),aspect.ratio=1)
    #########保存数据#########
    alldata$class <- str_replace_all(alldata$class, c("change-up" = paste0("protein-down&",ptmn,"-up"),
                                                      "more-up" = paste0(ptmn,"-more up"),
                                                      "change-down" = paste0("protein-up&",ptmn,"-down"),
                                                      "more-down" = paste0(ptmn,"-more down")))
    savexlsx(alldata,paste0(savepath,"scatterplot-",cg[i],".xlsx"),sheet = cg[i])
    ggplotsave(plot = pp,savepath=savepath,
               mapname = paste0("scatterplot-",cg[i]),
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily,
               ...) 
    diffdata <- alldata[alldata$class != "not-significant",]
    prosample <- prodata[,c("Accession",intersect(samp$Sample,names(prodata)))]
    names(prosample)[-1] <- paste(names(prosample)[-1],"protein",sep = "_")
    site[,"prosite"]<-gsub(";.*","",site$`Positions within proteins`)%>%paste0(site$Accession,":",site$`Amino acid`,.)
    sitesample <- site[,c("prosite",intersect(samp$Sample,names(site)))]
    names(sitesample) <- c("prosite",paste(names(sitesample)[-1],paste0(ptmn,"protein"),sep = "_"))
    left_join(diffdata,sitesample,by="prosite") %>% left_join(.,prosample,by="Accession") %>% left_join(.,anno,by="Accession") -> diffalldata
    wb <- addsheet(data = diffalldata,wb = wb,sheet = cg[i])
    addStyle(wb, sheet = cg[i], style = createStyle(fgFill = "#FBEAFF"), rows=(which(diffalldata$class == paste0("protein-down&",ptmn,"-up"))+1),cols = which(names(diffalldata)=="class"))
    addStyle(wb, sheet = cg[i], style = createStyle(fgFill = "#f47d8c"), rows=(which(diffalldata$class == paste0(ptmn,"-more up"))+1),cols = which(names(diffalldata)=="class"))
    addStyle(wb, sheet = cg[i], style = createStyle(fgFill = "#7edbfc"), rows=(which(diffalldata$class == paste0("protein-up&",ptmn,"-down"))+1),cols = which(names(diffalldata)=="class"))
    addStyle(wb, sheet = cg[i], style = createStyle(fgFill = "#727fb5"), rows=(which(diffalldata$class == paste0(ptmn,"-more down"))+1),cols = which(names(diffalldata)=="class"))
    savewb(wb = wb,filename = paste0("去除本底差异数据.xlsx"),overwrite = T)
    
  }
  print("~联合分析散点图绘制完成~")
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 10, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 8, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-v","--vesion", default = "差异表达矩阵(未筛选)",help = "数据表名称，默认为差异表达矩阵(未筛选)")
  parser$add_argument("-pf","--proFC", default = NULL,type= "double",help = "蛋白组FC值，默认读取rawdata路径下的ratio文件")
  parser$add_argument("-sf","--phosFC", default = NULL,type= "double",help = "修饰组FC值，默认读取rawdata路径下的ratio文件")
  parser$add_argument("-p","--pvalue", default = NULL,type= "double",help = "p值，默认读取rawdata路径下的ratio文件")
  parser$add_argument("-ip","--inputpath", default = "rawdata/",help = "输入文件路径，默认为rawdata/")
  parser$add_argument("-sp","--savepath",default = "Scatterplot/", help = "输出结果路径，默认为当前路径")
  parser$add_argument("-bp","--backpath", default = "background/",help = "Annotation_Data文件路径，默认为.background/")
  args <- parser$parse_args()
  unionvol <- do.call(unionvolplot,args = args)
}
