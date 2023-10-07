#' transmetaptp
#' 空转空代点对点联合分析
#' 
#' @param species 物种
#' @param samplename 样本名称
#' @param moderange 离子模式

#' @export
transmetaptp<-function(
    species = NULL,
    samplename=NULL,
    moderange=c("neg","pos"),
    ...
    
){

  if(is.null(samplename)){
    data=readdata(filename = "项目登记单.xlsx", sheet = "样本信息")
    samplename=data$META
  }
  
  if(is.null(species)){
    data1<-readdata(filename = "项目登记单.xlsx", sheet = "项目登记单")
    species=data1[data1[, 1] == "映射简写", 2]
  }  
  for(sample in samplename){
    for(mode in moderange){
      imzmltotransdata(samplename = sample, mode = mode,...)
      getuniondata(samplename= sample , mode =mode,...)
      TransAndMetaboanalyst(samplename = sample, mode =mode,savapath = paste0("./sample/union/result/",mode,"/"),...)
    }
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-spe","--species",default =NULL,help = "物种")
  parser$add_argument("-sn","--samplename",default =NULL,help = "样本名称")
  parser$add_argument("-mr","--moderange",default = c("neg","pos"), help = "离子模式")
  
  args <- parser$parse_args()
  
  
  result <- do.call(what = transmetaptp,args = args)
  
  
  
}
