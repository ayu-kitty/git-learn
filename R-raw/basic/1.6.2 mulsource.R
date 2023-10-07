#!/opt/conda/bin/Rscript

#' 加载包
#'
#' @param filepath 包路径
#' @param ... 见[source()]
#'
#' @export
mulsource <- function(filepath = "./",
                      ...) {
  
  file <- list.files(path = filepath, pattern = "\\.R$", full.names = T)

  for (i in file) {
    source(file = i, ...)
  }
}


#' 安装包
#'
#' @param filepath 包路径
#'
#' @export
roxygenizetoinstall <- function(filepath = "/home/lujw/lmbio",
                                build_vignettes = F) {
  rm(list = ls(pos = 1,all.names = T),pos = 1)
  
  rawpath <- getwd()
  setwd(filepath)
  filepath <- getwd()

  path <- list(package = filepath,
               database = "/data/hstore1")
  lmbioversion <- "2023"
  
  if(!dir.exists(paste0(filepath,"/data/"))){
    dir.create(paste0(filepath,"/data/"), recursive = T, mode = "0777")
    Sys.chmod(paths = paste0(filepath,"/data/"),mode = "0777",use_umask = F)
  }
  
  save(path,lmbioversion,file = paste0(filepath,"/data/path.rdata"))
  
  setwd("vignettes")
  
  thisfile <- list.files()
  
  for ( j in 1:length(thisfile)) {
    if(!file.exists(thisfile[j])){
      unlink(thisfile[j])
    }
  }
  
  dirpath <- list.dirs(path = "..",full.names = T,recursive = F)
  dirpath <- dirpath[!grepl(pattern = "vignettes$",x = dirpath)]
  dirpath <- dirpath[!grepl(pattern = "\\/\\.",x = dirpath)]
  allfile <- list.files(pattern = "\\.Rmd$",path =  dirpath,full.names = T,recursive = T)
  allfile <- allfile[!fs::link_exists(allfile)]
  allfilebase <- basename(allfile)
  allfilebase2 <- allfilebase[duplicated(allfilebase)]
  if(length(allfilebase2) > 0){
    stop(paste(allfilebase2,"文件重复",collapse = ";"))
  }
  
  if(length(allfile) > 0){
    for ( j in 1:length(allfile)) {
      
      if(fs::link_exists(allfilebase[j])){
        unlink(x = allfilebase[j]) 
      }
      
      if(!file.exists(allfilebase[j])){
        file.symlink(from = allfile[j],to = allfilebase[j])
      }else{
        stop(paste0(allfilebase[j],"在vignettes中有重复文件"))
      }
    }
  }
  
  setwd("../")
  
  setwd("R")

  thisfile <- list.files()
  
  for ( j in 1:length(thisfile)) {
    if(!file.exists(thisfile[j])){
      unlink(thisfile[j])
    }
  }
  
  allfile <- list.files(pattern = "\\.R$",path = "../R-raw",full.names = T,recursive = T)
  allfile <- allfile[!fs::link_exists(allfile)]
  
  allfilebase <- basename(allfile)
  allfilebase2 <- allfilebase[duplicated(allfilebase)]
  if(length(allfilebase2) > 0){
    stop(paste(allfilebase2,"文件重复",collapse = ";"))
  }
  
  if(length(allfile) > 0){
    for ( j in 1:length(allfile)) {

      if(fs::link_exists(allfilebase[j])){
        unlink(x = allfilebase[j]) 
      }

      if(!file.exists(allfilebase[j])){
        file.symlink(from = allfile[j],to = allfilebase[j])
      }else{
        stop(paste0(allfilebase[j],"在Rlink中有重复文件"))
      }
    }
  }
    
  setwd("../")
  
  setwd("command")
  
  thisfile <- list.files()
  
  for ( j in 1:length(thisfile)) {
    if(!file.exists(thisfile[j])){
      stop(paste0("command中",thisfile[j],"链接缺失"))
    }
  }
  
  setwd("../")
  
  roxygen2::roxygenize(package.dir = filepath)
  
  setwd("R")
  
  allfile <- list.files(pattern = "\\.R$",path = "./",full.names = T,recursive = T)

  funname <- NULL

  for ( i in 1:length(allfile)) {

    source(allfile[i])

    funname2 <- ls(pos = 1,all.names = T)
    funname2 <- funname2[funname2 != "allfile" &
                           funname2 != "funname" &
                           funname2 != "i" &
                           funname2 != "funname2" &
                           funname2 != "roxygenizetoinstall" &
                           funname2 != "rawpath"]
    
    if(any(funname2 %in% funname)){
      dupname <- funname2[funname2 %in% funname]
      stop(paste0(allfile[i],"文件中",paste(dupname,collapse = ";"),"函数名重复"))
    }else{
      funname <- c(funname,funname2)
      rm(list = funname2,pos = 1)
    }
  }
  
  # funname3 <<- funname
  setwd("../")
  
  devtools::install(pkg = filepath,
                    build_vignettes = build_vignettes,
                    upgrade = "never",force = T) 

  setwd(rawpath)
}

#' gitlab安装
#'
#' @export
gitlabinstall <- function(){
  devtools::install_git(url = "http://lujiawei:oekey-gpQKemKzKzw7f1c28gsM@gitlab.oebiotech.cn/lm/lmbio.git",
                        build_vignettes = T,upgrade = "never",force = T)
}

