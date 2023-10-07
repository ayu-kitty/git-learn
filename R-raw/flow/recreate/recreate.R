#!/opt/conda/bin/Rscript

#' recreate
#'
#' 售后重出报告
#'
#' @param project project项目编号
#' @param path path 获取数据路径
#' @param copydata copydata 是否复制原始报告中的数据矩阵
#' @param renamemap renamemap 是否修改色谱图名称/删除对应样本的色谱图
#' @export
recreate <- function(project = "-aa$",
                     S3 = "OutData",
                     path = "../",
                     copydata = T,
                     renamemap = T,
                     cxlsx = F,
                     fram = "数据矩阵.xlsx",
                     fram_lc = "数据矩阵-LC.xlsx",
                     fram_gc = "数据矩阵-GC.xlsx",
                     ...) {
  UseMethod("recreate")
}


recreate.default <- function(project = "-aa$",
                             S3 = "OutData",
                             path = "../",
                             copydata = T,
                             renamemap = T,
                             cxlsx = F,
                             fram = "数据矩阵.xlsx",
                             fram_lc = "数据矩阵-LC.xlsx",
                             fram_gc = "数据矩阵-GC.xlsx",
                             ...){
  suppressMessages(library("plotly"))
  ###复制内部分析单
  if(cxlsx){
    file.copy("/data/hstore4/database/mould/rename/rename.xlsx",to="./")
  }else{
    name<-list.files(path = path,pattern=project)
    file<-list.files(path =  paste0(path,name,"/"),pattern=c("内部分析单*"),full.names = F,recursive = T)
    file.copy(from = paste0(path,name,"/",file) ,to="./")
    
    dataall<-lmbio::readdata(filename="内部分析单.xlsx",sheet="分析基本信息")
    rownames(dataall)<-dataall$key
    if(dataall["S3",2] %in% c("UntargetBothEmdb","UntargetBothPmdb","UntargetBoth")){
      #双平台
      #da_all<-lmbio::readdata(filename="内部分析单.xlsx",sheet="分析基本信息")
      #rownames(da_all)<-da_all$key
      #da_all["S3",2]=S3
      #da_all["数据矩阵",2]=fram
      #lmbio::savexlsx1(da_all,"内部分析单.xlsx",sheet="分析基本信息")
      ##lc
      data_lc<-lmbio::readdata(filename="内部分析单-lc.xlsx",sheet="分析基本信息")
      rownames(data_lc)<-data_lc$key
      data_lc["S3",2]=S3
      data_lc["数据矩阵",2]=fram_lc
      data_lc["missvalue",2]="NA"
      data_lc["rsd",2]="NA"
      data_lc["zeroprocess",2]="NA"
      lmbio::savexlsx1(data_lc,"内部分析单-lc.xlsx",sheet="分析基本信息")
      ##gc
      data_gc<-lmbio::readdata(filename="内部分析单-gc.xlsx",sheet="分析基本信息")
      rownames(data_gc)<-data_gc$key
      data_gc["S3",2]=S3
      data_gc["数据矩阵",2]=fram_gc
      data_gc["missvalue",2]="NA"
      data_gc["zeroprocess",2]="NA"
      lmbio::savexlsx1(data_gc,"内部分析单-gc.xlsx",sheet="分析基本信息")
    }else{
      data<-lmbio::readdata(filename="内部分析单.xlsx",sheet="分析基本信息")
      rownames(data)<-data$key
      data["S3",2]=S3
      data["数据矩阵",2]=fram
      data["missvalue",2]="NA"
      data["zeroprocess",2]="NA"
      if(!is.na(data["rsd",2])){
        data["rsd",2]="NA"
      }
      lmbio::savexlsx1(data,"内部分析单.xlsx",sheet="分析基本信息")
    }
    ###复制数据矩阵
    if(copydata){
      sj<-list.files(dir(paste0(path,name,"/",dataall["项目报告",2]),pattern = "数据矩阵$",full.names = T),full.names = T)
      file.copy(from = sj ,to="./")
    }
    if(file.exists("数据矩阵.xlsx")){
      if(dataall["S3",2] %in% c("UntargetBothEmdb","UntargetBothPmdb","UntargetBoth")){
        file.rename(from = "数据矩阵.xlsx",to = "数据矩阵-old.xlsx")
        dataframe <-lmbio::readdata(filename = "数据矩阵-old.xlsx",sheet="数据矩阵") 
        dataframe <- select(dataframe,-c("Annotation","ID Annotation"))
        lmbio::savexlsx1(dataframe,"数据矩阵.xlsx",sheet="数据矩阵")
        dataframe_lc <-dataframe[dataframe$`Ion mode` %in% c("LCMS-neg","LCMS-pos"),]
        dataframe_lc$`Ion mode` <- gsub("LCMS-pos","pos",dataframe_lc$`Ion mode`)
        dataframe_lc$`Ion mode` <- gsub("LCMS-neg","neg",dataframe_lc$`Ion mode`)
        lmbio::savexlsx1(dataframe_lc,"数据矩阵-LC.xlsx",sheet="数据矩阵")
        dataframe_gc<-dataframe[dataframe$`Ion mode` %in% c("GCMS"),]
        lmbio::savexlsx1(dataframe_gc,"数据矩阵-GC.xlsx",sheet="数据矩阵")
        
      }else{
        file.rename(from = "数据矩阵.xlsx",to = "数据矩阵-old.xlsx")
        dataframe <-lmbio::readdata(filename = "数据矩阵-old.xlsx",sheet="数据矩阵")
        dataframe <- select(dataframe,-c("Annotation","ID Annotation"))
        lmbio::savexlsx1(dataframe,"数据矩阵.xlsx",sheet="数据矩阵")
      }
    }
    ###链接raw
    if (!dir.exists("./raw/")) {
      if(dir.exists(paste0(path,name,"/raw/"))){
        system(command = paste0("ln -s ",path,name,"/raw/ ", "raw"))
      }else if(dir.exists(paste0(path,"raw/"))){
        system(command = paste0("ln -s ",path,"raw/ ", "raw"))
      } else{
        print("~raw目录不存在，请确认")
      }
    }
    ###复制/生成init.txt  
    if(file.exists(paste0(path,name,"/init.txt"))){
      file.copy(from = paste0(path,name,"/init.txt") ,to="./")
    }else{
      file.copy(from = paste0(path,name,"/分析确认单.xlsx") ,to="./")
      lmbio::ArrangeInfo(overwrite = F)}
    
    ###复制色谱图
    lmbio::createdir(filename = "2.色谱图",force = force)
    if(renamemap){    
      rename<-lmbio::readdata(filename="rename.xlsx",sheet=1)
      if("GCMS" %in% list.files(path= paste0(path,name,"/",dataall["项目报告",2],"/2.色谱图"))){
        remap_gclc(oldpath_lc=paste0(path,name,"/",dataall["项目报告",2],"/2.色谱图/LCMS/"),
                   oldpath_gc=paste0(path,name,"/",dataall["项目报告",2],"/2.色谱图/GCMS/"),
                   path_lc="./2.色谱图/LCMS/",
                   path_gc="./2.色谱图/GCMS/",
                   rename = rename)
      }else{
        remap(oldpath=paste0(path,name,"/",dataall["项目报告",2],"/2.色谱图/"),
              newpath="./2.色谱图/",rename = rename)
      }
    }else{
      move_map(mappath=paste0(path,name,"/",dataall["项目报告",2]),
               movepath="./2.色谱图/")
    }
    
    ###链接background
    if(!file.exists("./background")){
      system(command = paste0("ln -s ",path,name,"/background/ ", "background"))
    }
    ###链接实验报告文件夹
    filename<-paste0(path,name,"/实验报告/")
    system(command = paste0("ln -s '",filename,"' '","./","'"))
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-pro","--project",default = "-aa$", help = "项目编号")
  parser$add_argument("-p","--path",default = "../", help = "获取数据路径")
  parser$add_argument("-cx","--cxlsx", default = F,  action = "store_true", help="是否创建rename表格，默认不创建")
  parser$add_argument("-cd","--copydata", default = T,  action = "store_false", help="是否复制报告中的数据矩阵,加上参数不复制")
  parser$add_argument("-r","--renamemap", default = T,  action = "store_false", help="是否修改色谱图名称/删除色谱图，如果参数打开需准备rename.xlsx文件用于改图片名称")
  parser$add_argument("-f","--fram",default = "数据矩阵.xlsx", help = "数据矩阵")
  parser$add_argument("-fl","--fram_lc",default = "数据矩阵-LC.xlsx", help = "数据矩阵-LC")
  parser$add_argument("-fg","--fram_gc",default = "数据矩阵-GC.xlsx", help = "数据矩阵-GC")
  parser$add_argument("-S","--S3",default = "OutData", help = "项目类型")
  args <- parser$parse_args()
  
  result <- do.call(what = recreate,args = args) 
  
}
