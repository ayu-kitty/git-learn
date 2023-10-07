#! /opt/conda/bin/Rscript

#' @export
flow_autodock_Analyze <- function(name = "HXD-疾病交集靶点.xlsx",
                                  project = "",
                                  result = ifelse(project == "","分子对接分析",paste(project,"分子对接分析",sep = "-")),
                                  smiles = "Ingredient_Smiles",
                                  protein = "UniProtID",
                                  maxnum = NULL,
                                  mode = "LM",
                                  CLOUD = F){
  wd <- getwd()
  
  data <- readdata(filename = name,filetype = c("xlsx","xls","txt","csv"))
  
  if(!is.null(maxnum)){
    if(dim(data)[1] > maxnum){
      print(paste0("~最大运行数量为",maxnum))
      data <- data[1:maxnum,]
    }
  }
  
  data <- data[,c(smiles,protein),drop = F]
  data <- data[!is.na(data[,smiles]),,drop = F]
  data <- data[data[,smiles] != "",]
  data <- data[!is.na(data[,protein]),,drop = F]
  data <- data[data[,protein] != "",]
  data[,"resultfilename"] <- as.character(1:dim(data)[1])
  
  setwddir(filename = result)
  
  for ( i in 1:dim(data)[1]) {
    setwddir(data$resultfilename[i])
    system(command = paste0("source /etc/profile;source ~/.bashrc;",
                            packagepath(path = "command/flow_autodock_single"),
                            " '",data[i,smiles],"' ",
                            " '",data[i,protein],"' "))
    
    if(file.exists("ligand_out.pdbqt")){
      result <- readLines("ligand_out.pdbqt")
      result <- result[grepl(pattern = "^REMARK VINA RESULT:",x = result)]
      result <- gsub(pattern = "^REMARK VINA RESULT: *",replacement = "",x = result)
      result <- gsub(pattern = " +.+",replacement = "",x = result)
      result <- as.numeric(result)
	  #result <- as.numeric(unique(result)[1])#选择affinity 最低的值
      data[i,"affinity(kcal/mol)"] <- result
    }
    
    setwd("../")
  }
  
  data <- cbind(data[,c("resultfilename","affinity(kcal/mol)")],
                data[,!(colnames(data) %in% c("resultfilename","affinity(kcal/mol)"))])
  
  savexlsx1(data = data,filename = "分子对接结果.xlsx",sheet = "结果")
  
  # 生成报告
  mkreport(type = "分子对接",mode = mode,CLOUD = CLOUD)
  
  setwd(wd)
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-n","--name",default = "HXD-疾病交集靶点.xlsx", help = "文件名")
  parser$add_argument("-pr","--project",default = "", help = "项目编号-联系人")
  parser$add_argument("-r","--result",default = "分子对接分析", help = "结果目录")
  parser$add_argument("-s","--smiles",default = "Ingredient_Smiles", help = "smiles列名")
  parser$add_argument("-p","--protein",default = "UniProtID", help = "蛋白ID列名")
  parser$add_argument("-z","--zip",default = F, help = "是否压缩",action='store_true')
  parser$add_argument("-mn","--maxnum",default = NULL, type= "integer",help = "最大运行数量")
  parser$add_argument("-m","--mode",default = "LM", help = "公司模板")
  parser$add_argument("-cl","--CLOUD",default = F, help = "是否云交付",action='store_true')
  
  args <- parser$parse_args()
  zip <- args$zip
  args$zip <- NULL
  
  result <- do.call(what = flow_autodock_Analyze,args = args)
  
  if(zip){
    zip::zip(zipfile = "项目报告.zip",files = args$result)
    unlink(list.files(pattern = "[^z][^i][^p]$",include.dirs = T,all.files = T), recursive = T)
  }
  
}
