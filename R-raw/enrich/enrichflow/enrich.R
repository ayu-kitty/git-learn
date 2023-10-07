#!/opt/conda/bin/Rscript
#' 背景文件复制
#' @param db 背景文件夹路径
#'
#' @param dbfrom 背景库路径
#' @param keggmapfrom kegg背景库路径
#' @param org 物种缩写
#'
#' @export
movebackground <- function(db,
                           dbfrom,
                           keggmapfrom="/data/hstore3/database/kegg/",
                           org,
                           inputfile="F"){
  if(!file.exists(db)){
    copydir(from = paste0(dbfrom,"/",org),to = db)
    if(dbfrom != keggmapfrom){
      copydir(from = paste0(keggmapfrom,"/",org),to = db)
    }
  }
  
  progmt(inputpath = db,inputfile=inputfile)
}

#' 富集骨架函数
#'
#' @param type 数据库类型G/K/I/R/W
#' @param inputfile 差异diff文件,-diff-
#' @param savepath 富集结果保存路径
#' @param org 物种缩写
#' @param db 背景文件夹路径
#' @param fontfamily 绘图字体
#' @param ... 
#' @export
enrichpip<-function(type="K",inputfile = "A-B-diff-protein.xls",savepath = "./enrich/",org="hsa",
                    db="../background/", inputpath="./",
                    keggdb="/data/hstore3/database/kegg/",
                    fontfamily="sans",...){
  pacman::p_load(dplyr,openxlsx,this.path,stringr,Biostrings)
  compare<-gsub("-diff-.*","",inputfile)
  diffile<-paste0(inputpath,"/",inputfile)
  createdir(savepath)
  if(ncol(readdata(diffile))==1){
    fc="False"
  }else fc="True"
  
  if(type=="G"){
    result <- goenrich(db=db,con=diffile,outd=savepath) #生成富集分析统计表
    if(result == "stop"){
      savetxt(data = "未富集到条目",filename = paste0(savepath,"/GO/",compare,"/说明.txt"))
      return()
    }
    enrichbarplot(incompare = compare,savepath=savepath, type = type,filt = "T",fontfamily=fontfamily)#top30柱状图
    bubble(incompare = compare, savepath=savepath,type = type,filt = "T",grid = "T",fontfamily=fontfamily)#top15气泡图
    gotopbar(incompare = compare,savepath=savepath,fontfamily=fontfamily)#棒棒糖图
    mainchord(incompare = compare, savepath=savepath,type = type,filt = "T",fontfamily=fontfamily,inputfile = inputfile,inputpath = inputpath)#和弦图
    if(fc=="True"){
      butterfly(incompare = compare, savepath=savepath,type = type,fontfamily=fontfamily)#蝴蝶图
      circleplot(incompare = compare,savepath=savepath, type = type,fontfamily=fontfamily,inputfile = inputfile,inputpath = inputpath)#环状图
    }
    print(paste0("~",compare,"数据库GO富集绘图结束"))
  }else if(type=="K"){
    result <- keggenrich(db=db,con=diffile,outd=savepath) #生成富集分析统计表
    if(result == "stop"){
      savetxt(data = "未富集到通路",filename = paste0(savepath,"/KEGG/",compare,"/说明.txt"))
      return()
    }
    enrichbarplot(incompare = compare,savepath=savepath,type = type,filt = "T",number=10,fontfamily=fontfamily)#top20柱状图
    bubble(incompare = compare, type = type,savepath=savepath,filt = "T",grid = "F",fontfamily=fontfamily)#气泡图
    mainchord(incompare = compare, type = type,savepath=savepath,filt = "T",fontfamily=fontfamily,inputfile = inputfile,inputpath = inputpath)#和弦图
    if(fc=="True"){
      up_down(incompare = compare,savepath=savepath,fontfamily=fontfamily)#kegg上下调对比图
      butterfly(incompare = compare, savepath=savepath,type = type,filt = "T",fontfamily=fontfamily)#蝴蝶图
      circleplot(incompare = compare,savepath=savepath, type = type,fontfamily=fontfamily,inputfile = inputfile,inputpath = inputpath)#环状图
      system(paste0("keggmap -if '",inputfile,"' -f True -o ",org," -ip ",inputpath," -d ",db," -kd '",keggdb,"' -s '",paste0(savepath,"KEGG_map/'")))
    }else{
      system(paste0("keggmap -if '",inputfile,"' -f False -c gold -o ",org," -ip ",inputpath," -d ",db," -kd '",keggdb,"' -s '",paste0(savepath,"KEGG_map/'")))
    }
    print(paste0("~",compare,"数据库KEGG富集绘图结束"))
    
  }else if(type=="W" | type=="R"| type=="I"){
    result <- enrich(db=db,enclass=type,con=diffile,outd=savepath)  #生成富集分析统计表
    if(result == "stop"){
      savetxt(data = "未富集到通路",filename = paste0(savepath,"/",getenrichtype(type)$filn,"/",compare,"/说明.txt"))
      return()
    }
    enrichbarplot(incompare = compare,savepath=savepath, type = type,filt = "T",number=10,fontfamily=fontfamily)#top20柱状图
    bubble(incompare = compare, type = type,savepath=savepath,filt = "T",grid="F",number=20,fontfamily=fontfamily)#气泡图
    mainchord(incompare = compare, type = type,savepath=savepath,filt = "T",fontfamily=fontfamily,inputfile = inputfile,inputpath = inputpath)#和弦图
    if(fc=="True"){
      butterfly(incompare = compare, type = type,savepath=savepath,filt = "T",fontfamily=fontfamily)#蝴蝶图
    }
    print(paste0("~",compare,"数据库",getenrichtype(type)$filn,"富集绘图结束"))
  }
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  parser$add_argument("-t","--type",type="character", default="K", help="enrich database class,such as G/K/W/R/I", metavar="character")
  parser$add_argument("-if","--inputfile",default = "A-B-diff-protein.xls", help = "输入数据文件，数据文件名必须包含-diff-")
  parser$add_argument("-ip","--inputpath",default = "./", help = "输入文件路径,默认./")
  parser$add_argument("-s","--savepath",default = "enrich/", help = "富集结果保存路径,默认./enrich/")
  parser$add_argument("-kd","--keggdb",default = "/data/hstore3/database/kegg/", help="keggmap图所有物种存放位置，默认为/data/hstore3/database/kegg/")
  parser$add_argument("-o","--org",default = "hsa", type= "character",help = "物种kegg缩写,默认hsa")
  parser$add_argument("-db","--db",default = "../background/",help = "背景文件存放位置,默认../background/")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  args <- parser$parse_args()
  
  enrich_pipeline <- do.call(enrichpip,args = args)
}
