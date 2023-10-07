#!/opt/conda/bin/Rscript

#' 参数筛选函数
#' @export
filtpro <- function(filename = "Protein quantitation.xlsx",
                    path = "./rawdata/",
                    batchmethod = "Combat",
                    ...){
  fil <- paste0(path,filename)
  prodata <- readdata(fil)
  
  if("Score Sequest HT: Sequest HT" %in% colnames(prodata)){
    prodata <- filter(prodata,`Score Sequest HT: Sequest HT`>0)
  }
  if("# Unique Peptides" %in% colnames(prodata)){
    prodata <- filter(prodata,`# Unique Peptides`>0)
  }
  if("Unique peptides" %in% colnames(prodata)){
    prodata <- filter(prodata,`Unique peptides`>0)
  }
  if("Localization prob" %in% colnames(prodata)){
    prodata <- filter(prodata,`Localization prob`>=0.75)
  }
  if("Delta score" %in% colnames(prodata)){
    prodata <- filter(prodata,`Delta score`>=8)
  }
  
  if(file.exists("批次信息.xlsx")){
    print("~进行批次校正")
    prodata <- Batch_correct (input = prodata,
                              saminfo = "批次信息.xlsx",
                              savepath = paste0(path,"/批次校正"),
                              method = batchmethod)
  }
  
  savexlsx1(prodata,paste0(path,"/数据矩阵.xlsx"),overwrite = T)
  
  wd <- getwd()
  setwd(path)
  organizedata(rawdatafile = "数据矩阵.xlsx",
               rawclassfile = c("Sample_Group.xlsx","样品"),
               rawcomparefile = c("Sample_Group.xlsx","比较组"),
               datafile = "./oecloud/rawdata/rawdatafile.txt",
               classfile = "../classfile.yaml",
               classtypefile = "../classtype.xlsx", 
               keylist = 1,
               ...)
  setwd(wd)
}
#' 蛋白质组学处理下机数据函数
#'
#' @param inputpath 下机数据存储路径
#' @param savepath 处理后保存路径
#' @param fetch 是否进行下机数据处理
#'
#' @export
profetch <- function(inputpath="./projectdata/",savepath="./rawdata/",fetch="T",
                     missvarremovepercent = -1,
                     missvarremovebygroup = "auto",
                     missvarfillmethod = "auto",
                     rowNorm = "auto",
                     transNorm = "auto",
                     scaleNorm = "auto",
                     ref = "auto",
                     filter = "auto",
                     remainnum = -1,
                     qcFilter = "auto", 
                     rsd = -1,
                     ...){
  pacman::p_load(dplyr,openxlsx,stringr)
  wd <- getwd()
  createdir(savepath)
  if(fetch=="T"){
    if(!is.na(dir(path=inputpath,pattern="^Label set*")[1])){
      Mstype = "PDM"
    }else if(file.exists(paste0(inputpath,"MSMS.xlsx"))){
      Mstype = "PD"
    }else if(file.exists(paste0(inputpath,"proteinGroups.txt")) & file.exists(paste0(inputpath,"peptides.txt"))){
      Mstype = "LF"
    }else if(file.exists(paste0(inputpath,"modificationSpecificPeptides.txt"))){
      if(file.exists(paste0(inputpath,"样品标记对照表.xlsx"))){
        Mstype = "PT"
      }else Mstype = "PL"
    }else if(file.exists(paste0(inputpath,"外来数据.xlsx"))){
      Mstype = "PU"
    }else {Mstype = "SP"}
    
    write(x = Mstype,file = paste0(savepath,"/protein_type.txt"))
  }else{
    Mstype<-readdata(paste0(savepath,"/protein_type.txt"),header=F)[1,1]
  }
  
  switch ( Mstype,
           "PD" = {
             if(fetch=="T"){
               pdfetch(savepath=savepath)
             }
             prosamgroup(inputpath = savepath)
             filtpro(filename="Protein quantitation.xlsx",path=savepath,...)
             # TMT/iTRAQ单组标记
           },
           "PDM" = {
             if(fetch=="T"){
               pdmultifetch(savepath=savepath)
             }
             prosamgroup(inputpath = savepath)
             filtpro(filename="Protein quantitation.xlsx",path=savepath,...)
             # 多组标记PD搜库
           },
           "LF" = {
             if(fetch=="T"){
               lbfetch(savepath=savepath)
             }
             prosamgroup(inputpath = savepath)
             filtpro(filename="Protein quantitation.xlsx",path=savepath,...)
             # LF
           },
           "PT" = {
             if(fetch=="T"){
               ptfetch(savepath=savepath)
             }
             prosamgroup(inputpath = savepath)
             ptmfil<-dir(savepath,pattern = "*Site*")
             filtpro(filename=ptmfil,path=savepath,...)
             #修饰标记--磷酸化
           },
           "PL" = {
             if(fetch=="T"){
               plfetch(savepath=savepath)
             }
             prosamgroup(inputpath = savepath)
             ptmfil<-dir(savepath,pattern = "*Site*")
             filtpro(filename=ptmfil,path=savepath,...)
             #非标修饰
           },
           "SP" = {
             if(fetch=="T"){
               spfetch(savepath=savepath)
             }
             prosamgroup(inputpath = savepath)
             filtpro(filename="Protein quantitation.xlsx",path=savepath,...)
             # DIA(SP)
           },
           "PU" = {
             prosamgroup(inputpath = savepath)
             filtpro(filename="外来数据.xlsx",path=savepath,...)
             # 外来数据
           }
  )
  
  setwd(savepath)
  predata <- predealdataforfile(missvarremovepercent = getpredealparams(omic = "P",MStype = Mstype,name = "missvarremovepercent",value = missvarremovepercent),
                                missvarremovebygroup = getpredealparams(omic = "P",MStype = Mstype,name = "missvarremovebygroup",value = missvarremovebygroup),
                                missvarfillmethod = getpredealparams(omic = "P",MStype = Mstype,name = "missvarfillmethod",value = missvarfillmethod),
                                rowNorm = getpredealparams(omic = "P",MStype = Mstype,name = "rowNorm",value = rowNorm),
                                transNorm = getpredealparams(omic = "P",MStype = Mstype,name = "transNorm",value = transNorm),
                                scaleNorm = getpredealparams(omic = "P",MStype = Mstype,name = "scaleNorm",value = scaleNorm),
                                ref = getpredealparams(omic = "P",MStype = Mstype,name = " ref",value =  ref),
                                filter = getpredealparams(omic = "P",MStype = Mstype,name = "filter",value = filter),
                                remainnum = getpredealparams(omic = "P",MStype = Mstype,name = "remainnum",value =remainnum),
                                qcFilter = getpredealparams(omic = "P",MStype = Mstype,name = "qcFilter",value = qcFilter), 
                                rsd = getpredealparams(omic = "P",MStype = Mstype,name = "rsd",value = rsd),
                                classfile = "../classfile.yaml")
  
  predata <- select(predata,c(1,names(stylefun_group(sample = T,classfile = "../classfile.yaml"))))
  colnames(predata)[1] <- "Accession"
  savexlsx1(data = predata,filename = "绘图数据.xlsx")
  setwd(wd)
  return(Mstype)
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-f","--fetch",default = "T", help = "是否做下机数据标准提取，默认为T，否则只运行参数筛选与标准化")
  parser$add_argument("-ip","--inputpath",type="character", default="./projectdata/", help="下机数据文件存放路径，默认当前路径", metavar="character")
  parser$add_argument("-sp","--savepath",default = "./rawdata/", help = "处理后数据保存路径,默认当前路径")
  
  # 预处理参数
  parser$add_argument("-mrp","--missvarremovepercent",default = -1, type= "double",help = "缺失值筛选比例,默认0.5")
  parser$add_argument("-mrg","--missvarremovebygroup",default = "auto", help = "是否按组进行缺失值比例计算,使用参数将不按组进行计算")
  parser$add_argument("-mvf","--missvarfillmethod",default = "auto", 
                      help = "缺失值填充方式,包含none,mean_half,halfmin,valuemin,min,mean,median,knn_var,knn_smp,ppca,bpca,svdlmpute",
                      choices = c("none","mean_half","halfmin","valuemin","min","mean","median","knn_var","knn_smp","ppca","bpca","svdlmpute"))
  parser$add_argument("-rn","--rowNorm",default = "auto", 
                      help = "归一化方式,包含NULL,SumNorm,MedianMeanNorm,MedianNorm,QuantileNorm,SamplePQN(ref需提供样本名),GroupPQN(ref需提供组名),CompNorm(ref需提供蛋白或者代谢名),SpecNorm(ref需提供样本名)",
                      choices = c("NULL","SumNorm","MedianMeanNorm","MedianNorm","QuantileNorm","SamplePQN","GroupPQN","CompNorm"))
  parser$add_argument("-tn","--transNorm",default = "auto", help = "转化方式,包含LogNorm,Log2minNorm,SrNorm,CrNorm",
                      choices=c("NULL","LogNorm","Log2minNorm","SrNorm","CrNorm"))
  parser$add_argument("-sn","--scaleNorm",default = "auto", help = "标准化方式,包含MeanCenter,AutoNorm,ParetoNorm,RangeNorm",
                      choices = c("NULL","MeanCenter","AutoNorm","ParetoNorm","RangeNorm"))
  parser$add_argument("-r","--ref",default = "auto", help = "归一化方式提供参考样本或者特征的参数")
  parser$add_argument("-fi","--filter",default = "auto", 
                      help = "数量筛选方式,默认为mean,包含rsd,nrsd,mean,sd,iqr,mad,median,none",
                      choices = c("rsd","nrsd","mean","sd","iqr","mad","median","none"))
  parser$add_argument("-re","--remainnum",default = -1, type= "integer",help = "特征筛选上限,默认为100000,输入0为自动筛选数量")
  parser$add_argument("-qf","--qcFilter",default = "T",help = "是否进行QC样本的rsd筛选")
  parser$add_argument("-rs","--rsd",default = -1,type= "integer", help = "QC样本rsd筛选范围,默认30")
  
  args <- parser$parse_args()
  profetch_result <- do.call(profetch,args = args)
}
