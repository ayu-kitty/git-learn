#!/opt/conda/bin/Rscript

#' KEGG map上色函数

#' @param inputdir 包含-diff-文件的路径
#' @param savepath 保存路径
#' @param org 物种kegg缩写
#' @param inputfile -diff-具体文件名
#' @param db 背景表格文件夹
#' @param color 颜色
#' @param keggdb keggmap图文件夹，不含物种缩写
#' @param dbname gene_anno-kegg.backgroud.xls
#' @param ... 
#'
#' @export
keggmap <- function(inputpath="./",savepath = "./enrich/KEGG_map/",org="hsa",inputfile = "A-B-diff-protein.xls",db="../background/",
                    color = "crimson,RoyalBlue3,gold",fc = "True",keggdb = "/data/hstore3/database/kegg/",
                    dbname = "gene_anno-kegg.backgroud.xls",...){
  pacman::p_load(dplyr,R.utils)
  path0 = getwd()
  inputpath = getAbsolutePath(inputpath)
  savepath = getAbsolutePath(savepath)
  keggdb = getAbsolutePath(keggdb)
  db = getAbsolutePath(db)
  createdir(savepath)
  
  setwd(savepath)
  if(!file.exists(dbname)){
    system(paste0("cp ",db,"/",dbname," ./"))
  }
  if(!is.na(dir(inputpath,pattern = gsub("-diff-data","-list",inputfile))[1])){
    system(paste0("cp ",inputpath,"/",gsub("-diff-data","-list",inputfile)," ./"))
  }
  if(!file.exists(inputfile)){
    degDf<-read.delim(paste0(inputpath,"/",inputfile), sep="\t", header=T, quote="",check.names=FALSE)
    bad<-read.delim(paste0(db,"/",dbname), sep="\t", header=F, quote="",check.names=FALSE)
    degDf[,1]<-gsub(":.*","",degDf[,1])
    bad[,"id"]<-gsub(":.*","",bad[,1])
    degDf1<-merge(degDf,bad,by="id",sort=F)
    if(ncol(degDf)==3){
      degDf<-select(degDf1,c(4,2,3))
    }else degDf<-as.data.frame(degDf1[,2])
    names(degDf)[1]<-"id"
    savexls(degDf,inputfile,quote=F)
  }
  if( is.na(list.dirs(paste0(keggdb,"/html_png"))[1]) & is.na(dir(keggdb,pattern = ".html")[1])){
    mapdb<-paste0(keggdb,"/",org)
  }else if(!is.na(list.dirs(paste0(keggdb,"/html_png"))[1])){
    mapdb<-paste0(keggdb,"/html_png")
  }else mapdb<-keggdb
  if(fc=="False"){
    system(paste0("path_nodiff '",inputfile,"' ",dbname," ",mapdb," ",color,">/dev/null 2>&1"))
  }else{
    system(paste0("path_diff '",inputfile,"' ",dbname," ",mapdb," ",color,">/dev/null 2>&1"))
  }
  setwd(path0)
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  parser$add_argument("-if","--inputfile",default = "A-B-diff-protein.xls", help = "输入数据文件，数据文件名必须包含-diff-")
  parser$add_argument("-ip","--inputpath",default = "./", help = "-diff-文件所在路径，默认当前路径")
  parser$add_argument("-s","--savepath",default = "enrich/KEGG_map/", help = "富集结果保存路径,默认./enrich/KEGG_map/")
  parser$add_argument("-o","--org",default = "hsa", type= "character",help = "物种kegg缩写,默认hsa")
  parser$add_argument("-kd","--keggdb",default = "/data/hstore3/database/kegg/", help="keggmap图所有物种存放位置，默认为/data/hstore3/database/kegg/")
  parser$add_argument("-d","--db",default = "../background/",help = "背景文件存放位置,默认../background/")
  parser$add_argument("-dn","--dbname",default = "gene_anno-kegg.backgroud.xls",help = "背景文件名称,默认gene_anno-kegg.backgroud.xls")
  parser$add_argument("-f","--fc",default = "True", type= "character",help = "keggmap图是否显示上下调信息(True/False)，默认为True")
  parser$add_argument("-c","--color",default = "crimson,RoyalBlue3,gold", help="有上下调输入：“上调颜色,下调颜色,共有颜色”，分隔符为”,“，没有上下调输入：“颜色”，默认为red,blue,orange")
  
  args <- parser$parse_args()
  keggmap_pipeline <- do.call(keggmap,args = args)
}