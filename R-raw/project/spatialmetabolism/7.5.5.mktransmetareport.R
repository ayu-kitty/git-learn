#' mktransmetareport
#' 创建空转空代点对点联合分析报告

#' @param samplepath sample路径
#' @param registfile 登记单文件名
#' 
#' @export
mktransmetareport<-function(samplepath="../sample/",
                            registfile=NULL
    
){
  
  if(is.null(registfile)){
    registfile=list.files(path="./",pattern="确认单")
  }
  
  data <- readdata(filename = "项目登记单.xlsx", sheet = "项目登记单")
  sampletable <- readdata(filename = "项目登记单.xlsx", sheet = "样本信息")
  report <- paste0(data[data[, 1] == "项目编号", 2], "-", 
                   data[data[, 1] == "客户名称", 2], "-",
                   data[data[, 1] == "联系人名称", 2], 
                   "-空转空代联合-项目报告")
  setwddir(filename = report, force = T)
  copydir(from = paste0(samplepath,"union/map/Intensity/"),to = "./1.map/Intensity/")
  copydir(from = paste0(samplepath,"union/map/vague/"),to = "./1.map/vague/")
  copydir(from = paste0(samplepath,"union/uniondata/"),to = "./2.uniondata/")
  copydir(from = paste0(samplepath,"union/result/"),to = "./3.result/")
  if(length(which(file.exists(paste0(list.files(path="./3.result",full.names = T),"/2.pathway_map"))==TRUE))){
    notes<-"【2.pathway_map】文件夹中KEGG通路图上标红色的点为匹配到代谢物或者基因，颜色不代表上下调"
    write.table(notes,file="./3.result/说明.txt",row.names = FALSE,col.names = FALSE,quote = F)
  }
  file.copy(from=paste0("../",registfile),to="./4.样品登记单")
  file.rename(list.files(path="./4.样品登记单",pattern = ".xlsx",full.names = T),paste0("./4.样品登记单/样品登记单.xlsx"))
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-sp","--samplepath",default = "../sample/",help = "sample文件路径")
  parser$add_argument("-rf","--registfile",default = NULL, help = "样品登记单")
  
  args <- parser$parse_args()
  

  result <- do.call(what = mktransmetareport,args = args)


  
}