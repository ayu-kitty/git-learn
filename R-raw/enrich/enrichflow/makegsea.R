#!/opt/conda/bin/Rscript

#' 根据差异表格生成用于富集的文件
#' 
#' @param inputfile 差异文件名，默认差异代谢物.xlsx
#' @param inputpath 差异代谢物表格路径以及输出文件路径，默认当前路径
#' @param omic 组学类型，P/M，默认为M
#' @export
makegsea <- function(inputfile = "差异表达矩阵(未筛选).xlsx",
                     inputpath = "./",
                     outputpath = inputpath,
                     comparename = "",
                     groupfile="oecloud/rawdata/classfile.yaml",
                     omic = "M"){
  
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
  
  classdata <- readdata(filename = groupfile)
  
  if(omic=="M" | omic=="ML"){
    
    for(i in 1:length(com)){
      meta <- readdata(paste0(inputpath,"/",inputfile),sheet=com[i])
      
      if("diff_filter" %in% class(meta)){
        meta <- get_diffana_data_2(data = meta,datafrom = "diffdata")
      }
      
      com2 <- com[i]
      com2 <- unlist(strsplit(x = com2,split = "-vs-"))
      classdata2 <- classdata[com2]
      classdata3 <- unlist(classdata2)
      classdata4 <- NULL
      for ( j in 1:length(classdata2)) {
        classdata4 <- c(classdata4,rep(names(classdata2)[j],length(classdata2[[j]])))
      }
      classdata4 <- gsub(pattern = " ",replacement = "-",x = classdata4)
      com2 <- gsub(pattern = " ",replacement = "-",x = com2)
      
      needcol <- colnames(meta)[colnames(meta) %in% classdata3]
      
      if("Metabolites" %in% colnames(meta)){
        needcol <- c("Metabolites",needcol)
      }else{
        stop("无Metabolites列")
      }
      
      meta_1 <- select(meta,all_of(needcol))
      names(meta_1)[1] <- "id"
      
      meta_1 <- meta_1[!is.na(meta_1$id),]
      if(any(is.na(meta_1))){
        meta_1[is.na(meta_1)] <- min(as.matrix(meta_1[,-1]),na.rm = T)
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
      
      meta_1 <- meta_1[!duplicated(meta_1$id),]
      gseacls <- paste0(length(needcol)," 2 1\n",
                        paste0("#",paste(com2,collapse = " "),"\n"),
                        paste0("",paste(classdata4,collapse = " "),""))
      
      savexls(filename = paste0(outputpath,"/",com[i],"-gsea-data.xls"),
              data = meta_1,
              quote = F)
      
      write(x = gseacls,file = paste0(outputpath,"/",com[i],"-gsea-data.cls"))
    }
    
  }else if(omic=="P" | omic=="PF"){
    
    for(i in 1:length(com)){
      
      pro <- readdata(paste0(inputpath,"/",inputfile),sheet=com[i])
      
      if("diff_filter" %in% class(pro)){
        pro <- get_diffana_data_2(data = pro,datafrom = "diffdata")
      }
      
      if("Gene Name"%in%names(pro)){
        pro[,"Gene Name"] <- gsub("  ","",pro[,"Gene Name"])
        pro[,"Gene Name"][is.na(pro[,"Gene Name"])] <- ""
        pro[,"id"]<-paste0(pro[,"Accession"],":",pro[,"Gene Name"])
      }else{
        pro[,"id"]<-pro[,"Accession"]
      }
      
      com2 <- com[i]
      com2 <- unlist(strsplit(x = com2,split = "-vs-"))
      classdata2 <- classdata[com2]
      classdata3 <- unlist(classdata2)
      classdata4 <- NULL
      for ( j in 1:length(classdata2)) {
        classdata4 <- c(classdata4,rep(names(classdata2)[j],length(classdata2[[j]])))
      }
      classdata4 <- gsub(pattern = " ",replacement = "-",x = classdata4)
      com2 <- gsub(pattern = " ",replacement = "-",x = com2)
      
      needcol <- colnames(pro)[colnames(pro) %in% classdata3]
      
      pro_1 <- select(pro,c("id",needcol))
      pro_1 <- pro_1[!is.na(pro_1$id),]
      if(any(is.na(pro_1))){
        pro_1[is.na(pro_1)] <- min(as.matrix(pro_1[,-1]),na.rm = T)
      }
      
      gseacls <- paste0(length(needcol)," 2 1\n",
                        paste0("#",paste(com2,collapse = " "),"\n"),
                        paste0("",paste(classdata4,collapse = " "),""))
      
      savexls(filename = paste0(outputpath,"/",com[i],"-gsea-data.xls"),
              data = pro_1,
              quote = F)
      write(x = gseacls,file = paste0(outputpath,"/",com[i],"-gsea-data.cls"))
    }
    
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-f","--inputfile",default = "差异表达矩阵(未筛选).xlsx", help = "差异文件名，默认差异代谢物.xlsx")
  parser$add_argument("-p","--inputpath",default = "./", help = "表格路径，默认当前路径")
  parser$add_argument("-o","--outputpath",default = "", help = "输出文件路径，默认当前路径")
  parser$add_argument("-c","--comparename",default = "", help = "比较组名称")
  parser$add_argument("-om","--omic",default = "M", help = "组学类型，P/M，默认为M")
  parser$add_argument("-g","--groupfile",default = "oecloud/rawdata/classfile.yaml", help = "样本名称文件")
  
  args <- parser$parse_args()
  
  if(args$outputpath == ""){args$outputpath <- args$inputpath}
  
  diffback <- do.call(makegsea,args = args)
  
}
