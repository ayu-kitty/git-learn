#!/opt/conda/bin/Rscript

#' 蛋白售后专用绘图函数
#'
#' @param hm 
#' @param volcano 
#' @param pca 
#' @param venn 
#' @param inputpath 
#' @export
psplot<-function(hm=F,hms=T,color = "heatmapcol",volcano=F,pca=F,venn=F,inputpath="./",showname=F,pfilter=0.05,fcfilter=1.5){
  #当前路径数据表格为Heatmap.xlsx
  if(hm){
    fn<-paste0(inputpath,"/Heatmap.xlsx")
    comm<-paste0("map_common_heatmap -f ",fn," -cg 1 -sr -sc -cr -co ",color)
    if(hms){
      system(paste0(comm," -s ./Heatmap/cluster_none"))
      system(paste0(comm," -s ./Heatmap/cluster_samples -cc"))
    }else{
      system(paste0("map_common_heatmap -f ",fn," -co ",color," -cg 1 -sc -sr -s ./Heatmap"))
    }
  }
  if(volcano){
    if(showname){
      comm<-paste0("map_common_volcano -f volcano.xlsx -ff ",fcfilter," -pf ",pfilter," -vf 0 -sn -nl name -i {png,pdf} -s ./Volcano/")
    }else comm<-paste0("map_common_volcano -f volcano.xlsx -ff ",fcfilter," -pf ",pfilter," -vf 0 -i {png,pdf} -s ./Volcano/")
    system(comm)
  }
  if(pca){
    print("uncomplete")
  }
  if(venn){
    if(!is.na(grep(".rds",list.files(path = inputpath))[1])){
      system(paste0("map_common_venn -f ",inputpath,"/* -s ./Venn"))
    }else{
      fn<-grep("表达矩阵.xlsx",dir(path = inputpath,pattern = "*.xlsx$",full.names = T),invert = T,value = T)
      
      for(com in fn){
        system(paste0("map_common_venn -f ",com," -s ./Venn/",gsub(".xlsx","",basename(com))))
        if(file.exists(paste0(inputpath,"/表达矩阵.xlsx"))){
          vr<-readxlsx(paste0("./Venn/",gsub(".xlsx","",basename(com)),"/Venn.xlsx"))
          zd<-readxlsx(paste0(inputpath,"/表达矩阵.xlsx"))
          names(zd)[1]<-"ID"
          newvr<-merge(vr,zd,by="ID",all.x=T,sort=F)
          savexlsx(newvr,paste0("./Venn/",gsub(".xlsx","",basename(com)),"/Venn.xlsx"),sheet = "venn")
        }
      }
      
    }
    
  }
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  parser$add_argument("-pf","--pfilter", default =  0.05,type = "double",
                      help = "火山图pvalue的筛选标准")
  parser$add_argument("-ff","--fcfilter", default =  1.5,type = "double",
                      help = "火山图fc的筛选标准")
  parser$add_argument("-sn","--showname", default =  F, action = "store_true",
                      help = "火山图是否显示名称,默认不显示")
  parser$add_argument("-ip","--inputpath", default="./", help="输入文件路径，默认当前路径")
  parser$add_argument("-hm","--hm", default = F, action = "store_true", 
                      help = "是否绘制热图，热图数据文件名为Heatmap.xlsx，默认F")
  parser$add_argument("-co","--color",default = "heatmapcol",help = "热图配色设置")
  parser$add_argument("-hms","--hms", default = T, action = "store_false", 
                      help = "是否绘制热图蛋白聚类，默认T")
  parser$add_argument("-pca","--pca", default= F, action = "store_true", 
                      help = "是否绘制PCA，默认F")
  parser$add_argument("-ve","--venn", default= F, action = "store_true", 
                      help = "是否绘制venn，用输入文件路径下的所有rds或xlsx文件，如需拼接数据，来源表格名称必须为表达矩阵.xlsx，默认F")
  parser$add_argument("-vl","--volcano", default= F, action = "store_true", 
                      help = "是否绘制volcano，火山图数据文件名为volcano.xlsx，默认F")
  args <- parser$parse_args()
  
  psp <- do.call(psplot,args = args)
  
}
