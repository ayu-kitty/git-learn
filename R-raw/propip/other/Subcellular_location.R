#!/opt/conda/bin/Rscript
#' uniprot数据库ID亚细胞定位脚本
#' @export
subcelloc<-function(savepath = "./Subcellular_location/",inputfile="差异蛋白筛选结果.xlsx",unifile="uniprot.xlsx",savefile = "Subcellular_location.xlsx",imagetype = c("pdf","png") ,height = 8,width = 12,dpi=300,fontfamily="sans",...){
  pacman::p_load(openxlsx,readxl,ggplot2,dplyr)
  createdir(savepath)
  #细胞类型
  cellclass<-c("Lysosome","Centriole","Cytoplasm","Golgi apparatus",
               "Rough endoplasmic reticulum","Smooth endoplasmic reticulum membrane",
               "Mitochondrion","Cell membrane","Ribosome","Membrane",
               "Nuclear pore","Nucleus membrane","Nucleus","Nucleolus",
               "Vacuole","Synaose","Peroxisome","chloroplas",
               "Cytoskeleton","Microsome","Spindle pole","Endosome",
               "Cell wall","Secreted")
  grecell<-function(dat){
    if(is.na(dat)){
      NA
    }else{
      which(sapply(cellclass, function(x)grep(x,dat))==1) %>% names %>% paste0(.,collapse=";")
    }
  }
  uniprot <- read_excel(unifile)%>%as.data.frame()
  uniprot$`Subcellular location [CC]`<-gsub("Note=.*","",uniprot$`Subcellular location [CC]`) %>% gsub("SUBCELLULAR LOCATION: ","",.)
  #信息整理
  sulo<-uniprot[,"Entry"] %>% as.data.frame()
  sulo[,"Subcellular_location"]<-sapply(uniprot$`Subcellular location [CC]`,grecell)
  names(sulo)[1]<-"Accession"
  #鉴定蛋白匹配注释信息
  if(inputfile=="all"){
    fils<-list()
    fils[[1]]<-readxlsx("表达矩阵.xlsx")
    names(fils)[1]<-"credit"
    shs<-getsheetname("差异表达矩阵.xlsx")
    for(sh in shs){
      fils[[sh]]<-readxlsx("差异表达矩阵.xlsx",sheet=sh)
    }
    savexlsx(fils,"蛋白列表.xlsx")
    inputfile="蛋白列表.xlsx"
  }
  sheetna<-getSheetNames(inputfile)
  wb <- createWorkbook(paste0(savepath,savefile))
  for(i in 1:length(sheetna)){
    dat<-read_excel(inputfile,sheet=sheetna[i])%>%as.data.frame()
    datcelo<-merge(dat,sulo,by="Accession",all.x=T)
    sl<-na.omit(datcelo$Subcellular_location)
    sclp<-strsplit(sl,";") %>% unlist %>% table %>% as.data.frame
    sclp<-sclp[order(sclp[,2]),]
    sclp[,1]<-as.character(sclp[,1])
    
    if(nrow(sclp)>10){
      sclpp<-rbind(c("Others",sum(sclp[1:(nrow(sclp)-9),2])),tail(sclp,9))
      sclpp <- as.data.frame(sclpp,row.names = c(1:length(sclpp[,1])))
      as.numeric(sclpp$Freq)->sclpp$Freq
      label_value <- paste(round(sclpp$Freq/sum(sclpp$Freq) * 100, 1), '%', sep = '')
      sclpp["sl"]<-paste(sclpp[,1],"(",sclpp[,2],",",label_value,")")
      sclpp<-sclpp[order(sclpp$Freq),]
      sclpp[,3]<-factor(sclpp[,3],levels = sclpp[,3])
      ps<-ggplot(sclpp, aes(x = "sub", y = Freq, fill = sl)) + 
        geom_bar(stat = 'identity', position = 'stack',width = 1,col="white",size=1) +
        theme_bw()+
        scale_fill_manual(values=rev(SelectColors("tableau10medium",n = nrow(sclpp))))+
        coord_polar(theta = 'y')+
        labs(x = '', y = '', title = '',fill="")+
        guides(fill = guide_legend(reverse = TRUE))+
        theme(axis.title= element_blank(),axis.text = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank())+
        theme(legend.text=element_text(size=14))
    }else{
      sclpp <- as.data.frame(sclp,row.names = c(1:length(sclp[,1])))
      as.numeric(sclpp$Freq)->sclpp$Freq
      label_value <- paste(round(sclpp$Freq/sum(sclpp$Freq) * 100, 1), '%', sep = '')
      sclpp["sl"]<-paste(sclpp[,1],"(",sclpp[,2],",",label_value,")")
      sclpp<-sclpp[order(sclpp$Freq),]
      sclpp[,3]<-factor(sclpp[,3],levels = sclpp[,3])
      ps<-ggplot(sclpp, aes(x = "sub", y = Freq, fill = sl)) + 
        geom_bar(stat = 'identity', position = 'stack',width = 1,col="white") +
        theme_bw()+
        scale_fill_manual(values=rev(SelectColors("tableau10medium",n = nrow(sclpp))))+
        coord_polar(theta = 'y')+
        labs(x = '', y = '', title = '',fill="")+
        guides(fill = guide_legend(reverse = TRUE))+
        theme(axis.title= element_blank(),axis.text = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank())+
        theme(legend.text=element_text(size=14))
    }
    ggplotsave(plot = ps,savepath=paste0(savepath,sheetna[i]),
               mapname = "Subcellular_location",
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily)

    names(sclpp)[1]<-"Subcellular_location"
    write.xlsx(sclpp[order(sclpp[,1]),],paste0(savepath,sheetna[i],"/Subcellular_location.xlsx"))
    
    addWorksheet(wb,sheetna[i] )
    writeData(wb, sheet = i, datcelo)
  }
  saveWorkbook(wb,paste0(savepath,savefile), overwrite = TRUE)

}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("png","pdf"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 12, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 9, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-u","--unifile",type="character", default="uniprot.xlsx", help="uniprot 下载的亚细胞定位文件名", metavar="character")
  parser$add_argument("-if","--inputfile",type="character", default="差异蛋白筛选结果.xlsx", help="需要注释的文件，默认差异蛋白筛选结果.xlsx", metavar="character")
  parser$add_argument("-s","--savepath",default = "./Subcellular_location/", help = "结果保存路径,默认./Subcellular_location/")
  parser$add_argument("-sf","--savefile",default = "Subcellular_location.xlsx", help = "结果文件名,默认Subcellular_location.xlsx")
  args <- parser$parse_args()
  subcellocation <- do.call(subcelloc,args = args)
}






