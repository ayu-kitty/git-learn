#!/opt/conda/bin/Rscript


#' 绘图售后
#'
#' @param inputpath 数据路径
#' @param hm hm是否绘制heatmap
#' @param hms hms是否聚类
#' @param mapdata mapdata绘图数据
#' @param compare compare是否按照比较组绘制heatmap（设置分组信息颜色,若设置参数需提供classtype.xlsx）
#' @param color color热图颜色
#' @param classtype classtype热图分组信息颜色表格
#' @param volcano volcano是否绘制火山图
#' @param showname showname火山图是否展示代谢物名称
#' @param name name火山图展示代谢物名称
#' @param pfilter pfilter火山图p值筛选标准
#' @param fcfilter fcfilter火山图FC值筛选标准
#' @param vipfilter vipfilter火山图vip值筛选标准
#' @param venn venn是否绘制venn
#' @param pca pca是否绘制PCA
#' @param pcaall 绘制百分之95置信区间
#' @param pcacol 绘制分组置信区间（绘制区间边框,无填充）
#' @param pca3d 绘制3dPCA
#' @export
metamap <- function(hm = F,
                    hms = T,
                    inputpath = "./",
                    mapdata = "heatmap.xlsx",#"数据矩阵.xlsx","volcano.xlsx"
                    compare = F,
                    color = "heatmapcol",
                    classtype = "classtype.xlsx",
                    volcano = F,
                    showname = F,
                    name = NULL,
                    pfilter = 0.05,
                    fcfilter = 0,
                    vipfilter = 1,
                    venn = F,
                    pca = F,
                    pcaall = F,
                    pcafill = F,
                    pcacol = F,
                    pca3d = F,
                    ...){
  options(warn=-1)
  if(hm){
    if(compare){
      fn<-paste0(inputpath,mapdata)
      comm<-paste0("map_common_heatmap -f ",fn," -sr -sc -cr -cg 1 -mn heatmap -cf classtype -co ",color)
      if(hms){
        #横聚类
        system(paste0(comm, " -s " ,paste0(inputpath,"Heatmap/cluster_none") ))
        #横纵聚类
        system(paste0(comm," -s " ,paste0(inputpath,"Heatmap/cluster_samples ","-cc")))
      }else{
        system(paste0("map_common_heatmap -f ",fn," -co ",color," -cg 1 " ," -sc "," -sr ", " -s ", paste0(inputpath,"Heatmap")))
      }
    }else{
      fn<-paste0(inputpath,mapdata)
      comm<-paste0("map_common_heatmap -f ",fn," -cg 1 -sr -sc -cr -mn heatmap -co ",color)
      if(hms){
        #横聚类
        system(paste0(comm, " -s " ,paste0(inputpath,"Heatmap/cluster_none") ))
        #横纵聚类
        system(paste0(comm," -s " ,paste0(inputpath,"Heatmap/cluster_samples ","-cc")))
      }else{
        system(paste0("map_common_heatmap -f ",fn," -co ",color," -cg 1 " ," -sc "," -sr ", " -s ", paste0(inputpath,"Heatmap")))
      }
    }
  }
  
  if(volcano){
    fn <- paste0(inputpath,mapdata)
    #fn<-list.files(path = inputpath,pattern = ".xlsx")
    if(showname){
      system(paste0("map_common_volcano -f ", fn ," -ff ",fcfilter," -pf ",pfilter," -vf ", vipfilter ," -sn ", "T",  " -nl " , name , " -s ", paste0(inputpath,"Volcano")))
    }else {
      system(paste0("map_common_volcano -f ", fn ," -ff ",fcfilter," -pf ",pfilter," -vf ", vipfilter , " -s", paste0(inputpath,"Volcano")))
    }
  }
  
  if(venn){
    if(!is.na(grep(".rds",list.files(path = inputpath))[1])){
      system(paste0("map_common_venn -f ",inputpath,"/*.rds -s ./Venn"))
    }else{
      fn<-grep(mapdata,dir(path = inputpath,pattern = "*.xlsx$",full.names = T),invert = T,value = T)
      for(com in fn){
        system(command = paste0("map_common_venn -f ",fn," -s ./Venn/",gsub(".xlsx","",basename(fn))," -mn Venn "))
        if(file.exists(paste0(inputpath,mapdata))){
          vr<-lmbio::readxlsx(paste0("./Venn/",gsub(".xlsx","",basename(fn)),"/Venn.xlsx"))
          zd<-lmbio::readxlsx(paste0(inputpath,mapdata))
          names(zd)[1]<-"ID"
          newvr<-merge(vr,zd,by="ID",all.x=T,sort=F)
          lmbio::savexlsx(newvr,paste0("./Venn/",gsub(".xlsx","",basename(fn)),"/Venn.xlsx"),sheet = "venn")
        }
      }
    }
  }
  
  if(pca){
    if(!is.na(grep(".rds",list.files(path = inputpath))[1])){
      if(pcaall){
        system(paste0("map_mulstatistics_scoremap -f ",inputpath,"*.rds" , " -s " , paste0(inputpath,"PCA1")))
      }
      if(pcafill){
        system(paste0("map_mulstatistics_scoremap2 -f ",inputpath,"*.rds" , " -s " , paste0(inputpath,"PCA2")))
      }
      if(pcacol){
        system(paste0("map_mulstatistics_scoremap2 -f ",inputpath,"*.rds"," -col ","T"," -fl ","F" , " -s " , paste0(inputpath,"PCA3"))) 
      }
    }else{
      if(pcaall){
        system(paste0("map_mulstatistics_scoremap -f ",inputpath,"*.xls" , " -s " , paste0(inputpath,"PCA1")))
      }
      if(pcafill){
        system(paste0("map_mulstatistics_scoremap2 -f ",inputpath,"*.xls" , " -s " , paste0(inputpath,"PCA2")))
      }
      if(pcacol){
        system(paste0("map_mulstatistics_scoremap2 -f ",inputpath,"*.xls"," -col ","T"," -fl ","F" , " -s " , paste0(inputpath,"PCA3"))) 
      }
    }}
  if(pca3d){
    if(!is.na(grep(".rds",list.files(path = inputpath))[1])){
      system(paste0("map_mulstatistics_3dscoremap -f ",inputpath,"*.rds" , " -s " , paste0(inputpath,"3DPCA")))
    }else{system(paste0("map_mulstatistics_3dscoremap -f ",inputpath,"*.xls" , " -s " , paste0(inputpath,"3DPCA")))}
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-in","--inputpath", default = "./",help = "输入文件路径,默认当前路径")
  parser$add_argument("-hm","--hm", default = F, action = "store_true", help = "是否绘制热图,热图数据文件名为Heatmap.xlsx,默认F")
  parser$add_argument("-hms","--hms", default = T, action = "store_false", help = "是否绘制热图聚类,默认T")
  parser$add_argument("-mp","--mapdata", default = "heatmap.xlsx", help = "绘图数据")
  parser$add_argument("-cm","--compare", default = F, action = "store_true", help = "是否按照比较组绘制heatmap,默认F")
  parser$add_argument("-co","--color",default = "heatmapcol",help = "热图配色设置")
  parser$add_argument("-ct","--classtype",default = "classtype.xlsx", help = "热图分组颜色模板")
  parser$add_argument("-vl","--volcano", default= F, action = "store_true", help = "是否绘制volcano,火山图数据文件名为volcano.xlsx,默认F")
  parser$add_argument("-sn","--showname", default =  F, action = "store_true",help = "火山图是否显示名称,默认不显示")
  parser$add_argument("-na","--name", default =  "name",help = "火山图显示名称")
  parser$add_argument("-pf","--pfilter", default =  0.05,type = "double",help = "火山图pvalue的筛选标准")
  parser$add_argument("-ff","--fcfilter", default =  0,type = "double",help = "火山图fc的筛选标准")  
  parser$add_argument("-vf","--vipfilter", default =  1,type = "double",help = "火山图vip的筛选标准")
  parser$add_argument("-ve","--venn", default= F, action = "store_true", help = "是否绘制venn,用输入文件路径下的所有rds或xlsx文件,默认F")
  parser$add_argument("-pc","--pca", default= F, action = "store_true", help = "是否绘制PCA,默认F")
  parser$add_argument("-pl","--pcaall", default= F, action = "store_true", help = "绘制置信区间,默认T")
  parser$add_argument("-pil","--pcafill", default= F, action = "store_true", help = "绘制分组置信区间,绘制填充颜色,默认F")
  parser$add_argument("-po","--pcacol", default= F, action = "store_true", help = "绘制分组置信区间,绘制区间边框,无填充,默认F")
  parser$add_argument("-pd","--pca3d", default= F, action = "store_true", help = "是否绘制3dpca")
  
  args <- parser$parse_args()
  
  result <- do.call(what = metamap,args = args) 
}