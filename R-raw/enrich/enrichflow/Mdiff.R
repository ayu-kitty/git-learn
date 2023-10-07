#!/opt/conda/bin/Rscript

#' 根据差异表格生成用于富集的文件
#' 
#' @param inputfile 差异文件名，默认差异代谢物.xlsx
#' @param inputpath 差异代谢物表格路径以及输出文件路径，默认当前路径
#' @param omic 组学类型，P/M，默认为M
#' @export
makediff <- function(inputfile = "差异表达矩阵.xlsx",
                     inputpath = "./",
                     outputpath = inputpath,
                     comparename = "",
                     omic = "M",
                     chebi = F){
  
  pacman::p_load(dplyr,readxl,openxlsx)
  
  com <- getsheetname(filename = paste0(inputpath,"/",inputfile))
  com <- setdiff(com,c(grep("原始",com,value = T),grep("可信",com,value = T)))
  if("" %in% com){
    if(comparename != ""){
      com <- comparename
    }else{
      com <- gsub("\\..*$","",basename(inputfile))
    }
  }
  
  if(omic=="M" | omic=="ML"){
    
    for(i in 1:length(com)){
      meta <- readdata(paste0(inputpath,"/",inputfile),sheet=com[i])
      
      if("diff_filter" %in% class(meta)){
        meta <- get_diffana_data_2(data = meta)
      }
      
      if(chebi){
        
        if("ChEBI" %in% colnames(meta)){
          keggname <- "ChEBI"
        }else{
          if("Metabolites" %in% colnames(meta)){
            print("~无ChEBI注释列，根据Metabolites列进行注释")
            meta <- getmetainfo(data = meta,idlist = "Metabolites",needlist = "ChEBI")
            keggname <- "ChEBI"
          }else{
            stop("无ChEBI注释列,需添加ChEBI或者Metabolites列")
          }
        }
        
      }else{
        
        if("kegg"%in% colnames(meta)){
          keggname <- "kegg"
        }else if("KEGG"%in% colnames(meta)){
          keggname <- "KEGG"
        }else{
          if("Metabolites" %in% colnames(meta)){
            print("~无KEGG注释列，根据Metabolites列进行注释")
            meta <- getmetainfo(data = meta,idlist = "Metabolites",needlist = "KEGG")
            keggname <- "KEGG"
          }else{
            stop("无KEGG注释列,需添加KEGG或者Metabolites列")
          }
        }
        
      }
      
      if("FoldChange"%in%names(meta)){
        meta_1 <- select(meta,all_of(c(keggname,"FoldChange")))%>%na.omit()
        if(dim(meta_1)[1] == 0){
          return()
        }
        meta_1[,"type"] <- sapply(1:nrow(meta_1), 
                                  function(x){if(meta_1$FoldChange[x]>1){"Up"}else{"Down"}})
      }else if("FC"%in%names(meta)){
        meta_1 <- select(meta,all_of(c(keggname,"FC")))%>%na.omit()
        if(dim(meta_1)[1] == 0){
          return()
        }
        meta_1[,"type"] <- sapply(1:nrow(meta_1), 
                                  function(x){if(meta_1$FC[x]>1){"Up"}else{"Down"}})
        names(meta_1)[2]<-"FoldChange"
      }else {
        meta_1 <- select(meta,c(keggname))%>%na.omit()
      }

      if(dim(meta_1)[1] == 0){
        return()
      }
      
      if(!is.numeric(meta_1[,1])){
        meta_2 <- meta_1[0,]
        for ( j in 1:dim(meta_1)[1]) {
          meta_1_1 <- unlist(strsplit(split = ";\n",x = meta_1[j,1]))
          meta_1_1 <- unlist(strsplit(split = "; ",x = meta_1_1))
          if(length(meta_1_1) == 0){
            next
          }
          for ( k in 1:length(meta_1_1)) {
            if(meta_1_1[k] == ""){
              next
            }
            meta_1_1_1 <- meta_1[j,,drop = F]
            meta_1_1_1[1,1] <- meta_1_1[k]
            meta_2 <- rbind(meta_2,meta_1_1_1)
            # break
          }
        }
        meta_1 <- meta_2
      }
      
      names(meta_1)[1] <- "id"
      
      if(chebi){
        meta_1[,1]<-paste0("ChEBI-",meta_1[,1])
        savexls(filename = paste0(outputpath,"/",com[i],"-diff-data_chebi.xls"),
                data = meta_1,
                quote = F)
      }else{
        savexls(filename = paste0(outputpath,"/",com[i],"-diff-data.xls"),
                data = meta_1,
                quote = F)
      }
      
    }
    
  }else if(omic=="P" | omic=="PF"){
    
    for(i in 1:length(com)){
      
      pro <- readdata(paste0(inputpath,"/",inputfile),sheet=com[i])
      
      if("diff_filter" %in% class(pro)){
        pro <- get_diffana_data_2(data = pro)
      }
      
      if("Gene Name"%in%names(pro)){
        pro[,"Gene Name"] <- gsub("  ","",pro[,"Gene Name"])
        pro[,"Gene Name"][is.na(pro[,"Gene Name"])] <- ""
        pro[,"id"]<-paste0(pro[,"Accession"],":",pro[,"Gene Name"])
      }else{
        pro[,"id"]<-pro[,"Accession"]
      }
      
      if("FoldChange"%in%names(pro)){
        pro_1<-select(pro,c("id","FoldChange"))%>%na.omit()
        pro_1[,"type"]<-sapply(1:nrow(pro_1), 
                               function(x){if(pro_1$FoldChange[x]>1){"Up"}else{"Down"}})
      }else if("FC"%in%names(pro)){
        pro_1<-select(pro,c("id","FC"))%>%na.omit()
        pro_1[,"type"]<-sapply(1:nrow(pro_1), 
                               function(x){if(pro_1$FC[x]>1){"Up"}else{"Down"}})
        names(pro_1)[2]<-"FoldChange"
      }else{
        pro_1=select(pro,c("id"))%>%na.omit()
      }
      savexls(filename = paste0(outputpath,"/",com[i],"-diff-data.xls"),
              data = pro_1,
              quote = F)
    }
    
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-f","--inputfile",default = "差异表达矩阵.xlsx", help = "差异文件名，默认差异代谢物.xlsx")
  parser$add_argument("-p","--inputpath",default = "./", help = "差异代谢物表格路径以及输出文件路径，默认当前路径")
  parser$add_argument("-o","--outputpath",default = "", help = "差异代谢物表格路径以及输出文件路径，默认当前路径")
  parser$add_argument("-c","--comparename",default = "", help = "比较组名称")
  parser$add_argument("-om","--omic",default = "M", help = "组学类型，P/M，默认为M")
  parser$add_argument("-cb","--chebi",default = F, help = "是否为chebi",action='store_true')
  args <- parser$parse_args()
  
  if(args$outputpath == ""){args$outputpath <- args$inputpath}
  
  diffback <- do.call(makediff,args = args)
  
}
