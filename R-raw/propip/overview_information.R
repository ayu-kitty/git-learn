#!/opt/conda/bin/Rscript

#' @export
over_table <- function(rawdatapath = "rawdata",
                       projectpath = "项目报告",
                       savepath = "report/result/1.Project_information/overview_information.xlsx",
                       ...){
  pacman::p_load(readxl,stringr,openxlsx,scales)
  
  # ##project_title
  title <- readdata(paste0(rawdatapath,"/project_title.txt"),header =F) %>% as.character()
  
  ##样本数量、分组
  classfile<- readdata(paste0('classfile.yaml'))
  sample<- data.frame(group = c(names(classfile)))
  sample$samples<- sapply(1:length(classfile),function(x){
    paste0(classfile[x][[1]],collapse = ',')
  })
  sample$number<- sapply(1:length(classfile),function(x){
    length(classfile[x][[1]])
  })
  
  # ##物种及样本
  init <- readdata(paste0(rawdatapath,"/init.txt"),header =F)
  colnames(init)<- c('目录','信息')
  
  ##蛋白\位点数量
  protein_num <- readdata(paste0(rawdatapath,"/Identified_number.xlsx"))
  
  ##总蛋白注释结果
  annotationfile <- list.files(path = projectpath,pattern = "annotation.xlsx",full.names = T,recursive = T)
  data <- readdata(annotationfile[1])
  col <- colnames(data)
  proteins <- length(data$Accession)
  anno <- data.frame()
  back_names <- c("GO_id","pathway","eggNOG_id","PFAM_Family","InterPro_id","WikiPathway_id","Reactome_id")
  names(back_names) <- c('GO_term',"pathway_description",'eggNOG_function_classification',"PFAM_description","InterPro_ENTRY_NAME","WikiPathway_description","Reactome_description")
  background <- back_names[which(back_names%in%col)]%>%gsub("_.*","",.)%>%gsub("pathway","KEGG",.)
  for (i in 1:length(background)) {
    if(names(background[i])%in%col){
      num<- na.omit(data[,names(background[i])])%>%length()
	  if(!num==0){
      percent<- paste0(round(num/proteins*100,2),"%")
      func<-sapply(1:nrow(data), function(x){
        strsplit(as.character(data[x,names(background[i])]) ,"\\|")
      })%>%unlist()
      table<- table(func)%>%as.data.frame()
      top5 <- head(table[order(table$Freq, decreasing = T), ], 5)[,1]%>%as.character()
      top<- paste0(rep(1:5),".",top5)%>%paste(.,collapse = "; ")
      back_ta<- data.frame(background[i],num,percent,top)
	  colnames(back_ta) <- c("Database", "Number of comments", "Annotation rate","Top 5 features")
    }else{
	back_ta<- data.frame(background[i],num,"0",NA)
	colnames(back_ta) <- c("Database", "Number of comments", "Annotation rate","Top 5 features")
	}
    anno <- rbind(anno,back_ta)
  }}
  
  ##筛选标准ratio
  if(length(getSheetNames(paste0(rawdatapath,'/Sample_Group.xlsx')))>1){
    config <- readdata(paste0(rawdatapath,'/oecloud/config/config_diff.yaml'))
    if(!length(list.files(path = projectpath,pattern = "volcano*",full.names = T,recursive = T))==0){
      ratio<- data.frame('FC'=config[["params"]][["diff_filter"]][["fcfilter"]],
                         'Pvalue'=config[["params"]][["diff_filter"]][["pfilter"]])
      colnames(ratio)<- c('FC','P-value')
    }else{
      ratio<- data.frame('FC'=config[["params"]][["diff_filter"]][["fcfilter"]])
      colnames(ratio)<- c('FC')
    }
    ##差异筛选结果
    sample_group <- readdata(paste0(rawdatapath,"/Sample_Group.xlsx"),sheet = "比较组")
    sample_group[,1]<- gsub("/","-vs-",sample_group[,1])
    sa <- sample_group[,1]
    dataR <- data.frame()
    for (i in 1:nrow(sample_group)) {
      ##差异数量部分
      difffilename <- list.files(path = paste0(rawdatapath,"/oecloud/diff/"),pattern = paste0("^diff-data-",sa[i],"_.*.xls$"),full.names = T)
      if(length(difffilename)==0){
        dataR[i,1] <- sample_group[i,1]
        dataR[i,2] <- 0
        dataR[i,3] <- 0
        dataR[i,4] <- 0
      }else{
        diff <- readdata( difffilename[1])
        dataR[i,1] <- sample_group[i,1]
        dataR[i,2] <- dim(diff)[1]
        dataR[i,3] <- length(which(diff$FoldChange>=1))
        dataR[i,4] <- length(which(diff$FoldChange<=1))
      }
      if(!length(list.files(path = projectpath,pattern = "volcano*",full.names = T,recursive = T))==0){
        dataR$`t.test`[i]<- ifelse(sample_group$配对[which(sa[i]==sample_group$比较组)]==TRUE,'paired','unpaired')
      }
    }
    colnames(dataR)[1:4]<- c('Compare','All','upregulated','downregulated')
  }else{
    ratio<- NA
    dataR<- NA
  }
  
  
  if(file.exists(paste0(rawdatapath,"/样品标记对照表.xlsx"))){
    table <- readdata(paste0(rawdatapath,"/样品标记对照表.xlsx"))
	  if(length(grep("^delet", table$样品编号))>=1){
       table <-table[-grep("^delet", table$样品编号),]
      }else{
       table <-table
      }
    wb <- list(title,init,sample,protein_num,anno,ratio,dataR,table)
    names(wb)<- c('project_title','init','sample_group','原始数据统计','protein_anno','筛选条件','diff_number','样品标记对照表')
  }else{
    wb <- list(title,init,sample,protein_num,anno,ratio,dataR)
    names(wb)<- c('project_title','init','sample_group','原始数据统计','protein_anno','筛选条件','diff_number')
  }
  write.xlsx(wb, savepath)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  
  parser$add_argument("-rp","--rawdatapath",default="rawdata",help="运行路径")
  parser$add_argument("-pp","--projectpath",default="项目报告",help="项目路径")
  parser$add_argument("-sp","--savepath",default="report/result/1.Project_information/overview_information.xlsx",help="保存路径")
  args <- parser$parse_args()
  over_table <- do.call(over_table,args = args)
}
