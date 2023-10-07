#!/opt/conda/bin/Rscript

#' @export 
getQIdata <- function(data,
                      datapath = c("/data/hstore4/lumingos/2023-LC",
                                   "/data/hstore4/项目数据/LCMS/2022",
                                   "/data/hstore4/lumingos/project",
                                   "/data/hstore4/lumingos/项目归档/project",
                                   "/data/hstore4/lumingos/项目归档/项目数据/LCMS/2022",
                                   "/data/hstore4/项目数据/双平台/2022"),
                      copypath = "./",
                      ...){
  # data <- readdata(filename)
  if("项目编号/订单号" %in% colnames(data)){
    project <- data$`项目编号/订单号`
  }else if("订单号" %in% colnames(data)){
    project <- data$`订单号`
  }
  
  for ( i in 1:length(project)) {
    projectpath <- NULL
    for ( j in 1:length(datapath)) {
      projectpath2 <- list.files(path = datapath[j],pattern = project[i],full.names = T,include.dirs = T)
      projectpath <- c(projectpath,projectpath2)
    }
    
    if(length(projectpath) > 0){
      for ( j in 1:length(projectpath)) {
        projectidfile <- list.files(path = projectpath[j],pattern = "ID\\.csv$",full.names = T,recursive = T)
        projectmfile <- list.files(path = projectpath[j],pattern = "M\\.csv$",full.names = T,recursive = T)
        projectmspfile <- list.files(path = projectpath[j],pattern = "\\.msp$",full.names = T,recursive = T)
        projectfile <- c(projectidfile,projectmfile,projectmspfile)
        if(length(projectfile) > 0){
          file.copy(from = projectfile,to = copypath) 
        }else{
          print(paste0("~未找到:",project[i])) 
        }
      }
    }else{
      print(paste0("~未找到:",project[i])) 
    }
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  tool_getQIdata <- map_autodraw$new(getQIdata)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "项目文件",required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  
  # 基本参数
  parser$add_argument("-d","--datapath", default = c("/data/hstore4/lumingos/2023-LC",
                                                     "/data/hstore4/项目数据/LCMS/2022",
                                                     "/data/hstore4/lumingos/project",
                                                     "/data/hstore4/lumingos/项目归档/project",
                                                     "/data/hstore4/lumingos/项目归档/项目数据/LCMS/2022",
                                                     "/data/hstore4/项目数据/双平台/2022"), 
                      help = "数据检索路径",
                      nargs="+")
  parser$add_argument("-c","--copypath", default = "./", help = "保存路径")
  
  args <- parser$parse_args()
  
  result <- do.call(what = tool_getQIdata,args = args) 
}

#' @export
tool_getQIdata <- map_autodraw$new(getQIdata)$draw
