#!/opt/conda/bin/Rscript

#' 原始数据绘图
#'
#' @param savepath 
#' @param inputpath 
#' @param fontfamily 
#' @param imagetype 
#' @param ... 
#' @export
unionplot <- function(savepath = "./2.Qualitative/Statistics/",
                      samgpath = "./2.Qualitative/Identified/",
                      annpath = "./2.Qualitative/Annotation/",
                      dbpath = "../background/",
                      crepath = "./3.Reliable_results/",
                      inputpath = "./",
                      classfile = "../classtype.xlsx",
                      ...){
  
  if(file.exists(paste0(dbpath,"/Annotation_Data.xlsx"))){
    annorateplot(inputpath = inputpath,
                 savepath=annpath,
                 dbpath=dbpath,
                 ...)
  }
  
  
  multimarks <- dir(inputpath,pattern = "Label set.*")
  
  if(length(multimarks)==0){
    if(file.exists(paste0(inputpath,"/Identified_number.xlsx"))){
      identifiedplot(savepath = samgpath,
                     inputpath = inputpath,...)
    }
    profile <- dir(inputpath,pattern = "Pro.*")[1]
    pepfile <- dir(inputpath,pattern = "Pep.*")[1]
    sitefile <- dir(inputpath,pattern = ".*Sites.*")[1]
    if(!is.na(profile)){
      seqcoverage(savepath = savepath,inputpath = inputpath,inputfile = profile,...)
      MWplot(savepath = savepath,inputpath = inputpath,inputfile = profile,...)
      pepcountplot(savepath = savepath,inputpath = inputpath,inputfile = profile,...)
      if(is.na(sitefile)){
        abundrankplot(savepath = savepath,inputpath = inputpath,inputfile = profile,...) #表达丰度散点图，需要Sample表
      }
    }
    if(!is.na(pepfile)){
      peplengthplot(savepath = savepath,inputpath = inputpath,inputfile = pepfile,...)
    }
    if(!is.na(sitefile)){
      if(sitefile=="Phospho (STY)Sites.xlsx"){
        phostyplot(savepath = savepath,inputpath = inputpath,inputfile = sitefile)
      }
      prosite(savepath = savepath,inputpath = inputpath,inputfile = sitefile,...)
      abundrankplot(savepath = savepath,inputpath = inputpath,inputfile = sitefile,...) #表达丰度散点图，需要Sample表
    }
    if(!is.na(pepfile) & !is.na(sitefile)){
      pepsite(savepath = savepath,inputpath = inputpath,inputfile = pepfile,...)
    }
  }else{
    for(label in multimarks){
      savepath1 <- paste0(savepath,label,"/")
      identifiedplot(label=label,savepath = samgpath,mapname=paste0("Identified-number-",label),inputpath = inputpath,...)
      abundrankplot(savepath = savepath1,inputpath=inputpath,inputfile = paste0(label,"/Protein quantitation.xlsx"),...)
      seqcoverage(savepath = savepath1,inputpath = paste0(inputpath,"/",label,"/"),...)
      MWplot(savepath = savepath1,inputpath = paste0(inputpath,"/",label,"/"),...)
      pepcountplot(savepath = savepath1,inputpath = paste0(inputpath,"/",label,"/"),...)
      peplengthplot(savepath = savepath1,inputpath = paste0(inputpath,"/",label,"/"),...)
    }
  }
  
  
  if(file.exists(paste0(inputpath,"Sample_Group.xlsx"))){
    profile<-dir(inputpath,pattern = "Pro*")[1]
    sitefile<-dir(inputpath,pattern = "*Sites*")[1]
    if(!is.na(sitefile)){
      sampleproplot(inputfile = sitefile,inputpath = inputpath,savepath = samgpath,classfile=classfile,...)
      groupproplot(inputfile = sitefile,inputpath = inputpath,savepath = samgpath,classfile=classfile,...)
    }else{
      sampleproplot(inputfile = profile,inputpath = inputpath,savepath = samgpath,classfile=classfile,...)
      groupproplot(inputfile = profile,inputpath = inputpath,savepath = samgpath,classfile=classfile,...)
    }
  }
  print("~原始数据统计图绘制完成！") 
  plotfile <- dir(inputpath,pattern = "*绘图数据*")[1]
  sampfile <- dir(inputpath,pattern = "*Sample_Group*")[1]
  if(!is.na(plotfile) & !is.na(sampfile)){
    proboxplot(savepath=paste0(crepath,"/Boxplot/"),inputpath = inputpath,classfile=classfile,...)   #可信蛋白箱线图
    prodensityplot(savepath=paste0(crepath,"/Densityplot/"),inputpath = inputpath,classfile=classfile,...)   #可信蛋白密度图
    procorrtree(savepath=paste0(crepath,"/Sampletreeplot/"),inputpath = inputpath,classfile=classfile,...)   #可信蛋白树状图
    prorsd(savepath=paste0(crepath,"/RSD/"),inputpath = inputpath,classfile=classfile,...)   #可信蛋白RSD图
    samplecorr(savepath=paste0(crepath,"/Samplecorrplot/"),inputpath = inputpath,classfile=classfile,...)   #可信蛋白样本相关性图
  }
  print("~可信蛋白绘制完成！") 
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  
  # 此图参数
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sh","--dbpath", default = "../background/",help = "背景文件存放位置,默认为../background/")
  parser$add_argument("-sp","--savepath",default = "./2.Qualitative/Statistics/", help = "定性信息输出结果路径，默认./2.Qualitative/Statistics/")
  parser$add_argument("-sg","--samgpath",default = "./2.Qualitative/Identified/", help = "鉴定情况输出结果路径，默认./2.Qualitative/Identified/")
  parser$add_argument("-ap","--annpath",default = "./2.Qualitative/Annotation/", help = "注释结果输出路径，默认./2.Qualitative/Annotation/")
  parser$add_argument("-cp","--crepath",default = "./3.Reliable_results/", help = "可信分布输出结果路径，默认./3.Reliable_results/")
  parser$add_argument("-cf","--classfile",default = "../classtype.xlsx", help = "模板文件，默认../classtype.xlsx")
  args <- parser$parse_args()
  rawplot <- do.call(unionplot,args = args)
}