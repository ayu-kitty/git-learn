#' organizetmdata
#' 整理空转空代点对点分析的数据
#' 
#' @param transpath 空转报告路径
#' @param metapath 空代sample路径
#' @param moderange 离子模式
#' @param featurefile filtered_feature_bc_matrix.h5文件
#' @param positionfile tissue_positions_list.csv文件
#' @param clustdirectory graphclust目录


#' @export
organizetmdata<-function(
    transpath=NULL,
    metapath ="../sample",
    moderange=c("neg","pos"),
    featurefile="filtered_feature_bc_matrix.h5",
    positionfile="tissue_positions_list.csv",
    clustdirectory="gene_expression_graphclust"
    
){
  source("h5ToCSV.R")
  if(is.null(transpath)){
    transpath=list.files(pattern="Report",full.names = T)
  }
  
  data <- readdata(filename = "项目登记单.xlsx", sheet = "样本信息")
  transname=data$TRANS
  metaname=data$META
  dir.create("sample/union/transdata",recursive = TRUE)
  dir.create("sample/vague/",recursive = TRUE)
  dir.create("sample/qualitative/",recursive = TRUE)
  copydir(from = paste0(metapath,"/qualitative/Qualitative.xlsx"),to = "./sample/qualitative/")
  
  for(sample in metaname){
    trans<-transname[which(metaname==sample)]
    copydir(from = paste0(transpath,"/1.SpaceRanger/",trans,"/",featurefile),to = paste0("./sample/union/transdata/",sample))
    copydir(from = paste0(transpath,"/1.SpaceRanger/",trans,"/spatial/",positionfile),to = paste0("./sample/union/transdata/",sample))
    copydir(from = paste0(transpath,"/1.SpaceRanger/",trans,"/analysis/clustering/",clustdirectory,"/clusters.csv"),to = paste0("./sample/union/transdata/",sample))
    h5ToCSV(path=paste0("./sample/union/transdata/",sample),samplename=sample,filename=featurefile)
    file.rename(paste0("./sample/union/transdata/",sample,"/clusters.csv"),paste0("./sample/union/transdata/",sample,"/clusters_infor.csv"))
    cluster_info<-readdata(filename=paste0("./sample/union/transdata/",sample,"/clusters_infor.csv"))
    cluster_info$`sampleid`<-sample
    cluster_info$`group`<-sample
    colnames(cluster_info)[which(colnames(cluster_info)=="Cluster")]<-"clusters"
    write.csv(cluster_info,file=paste0("./sample/union/transdata/",sample,"/clusters_infor.csv"),row.names=FALSE)
    
    
    for(mode in moderange){
      copydir(from = paste0(metapath,"/final/",sample,"-",mode,".imzML"),to = "./sample/vague/")
      copydir(from = paste0(metapath,"/final/",sample,"-",mode,".ibd"),to = "./sample/vague/")
      transcrood<-readdata(filename=paste0("./sample/union/transdata/",sample,"/",positionfile),header = FALSE)
      transcrood<-transcrood[,c(1,4,3)]
      colnames(transcrood)<-c("Barcode","transx","transy")
      centerdata<-readdata(filename=paste0("center-",sample,".csv"))
      croodtable<-merge(centerdata,transcrood,by="Barcode")
      croodtable$transx<-as.numeric(croodtable$transx)
      croodtable$transy<-as.numeric(croodtable$transy)
      croodtable$metax<-NA
      croodtable$metay<-NA
      savexlsx1(croodtable,filename="crood.xlsx",sheet=paste0(sample,"-",mode))
    }
    
    copydir(from="crood.xlsx",to="./sample/union/")
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-tp","--transpath",default =NULL,help = "空转报告路径")
  parser$add_argument("-mp","--metapath",default ="../sample",help = "空代sample路径")
  parser$add_argument("-mr","--moderange",default = c("neg","pos"), help = "离子模式")
  parser$add_argument("-ff","--featurefile",default ="filtered_feature_bc_matrix.h5",help = "空转filtered_feature_bc_matrix.h5文件名称")
  parser$add_argument("-pf","--positionfile",default ="tissue_positions_list.csv",help = "空转tissue_positions_list.csv文件名称")
  parser$add_argument("-cd","--clustdirectory",default ="gene_expression_graphclust",help = "空转graphclust目录名称")
  
  args <- parser$parse_args()
  
  
  result <- do.call(what = organizetmdata,args = args)
  
  
  
}
