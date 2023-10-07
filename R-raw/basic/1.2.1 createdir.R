#!/opt/conda/bin/Rscript

#' 生成目录
#' 
#' @param filename 向量，文件夹名
#' @param force 逻辑，文件夹存在时是否强制删除重新创建，默认F
#'
#' @export 
createdir <- function(filename = "./",
                      linkdir = F,
                      linkinfo = "runinfo.txt",
                      force = F){
  tryCatch({
    
    if(length(filename) == 0){
      return("./")
    }
    
    if(filename %in% c("./",".")){
      return("./")
    }
    
    if(linkdir){
      filename <- gsub("\\/$","",filename)
      
      write(x = paste0("指向生成时间:",format(Sys.time(), "%Y%m%d")),
            file = linkinfo,
            append = T)
      
      if (file.exists(filename)) {
        print(paste0(filename, "文件夹原始链接重新创建链接"))
        
        if(!fs::link_exists(filename)){
          
          rawfile <- paste0(filename,"-",format(Sys.time(), "%Y%m%d"),"-raw")
          
          if(file.exists(rawfile)){
            rawfile <- paste0(rawfile,"-",RandomCode(n = 5))
          }
          
          moverawfile(rawpath = filename,newpath = rawfile)
          
        }else{
          rawfile <- fs::link_path(filename)
          rawfile <- paste0(dirname(filename),"/",rawfile)
          
          unlink(filename,
                 recursive = T,
                 force = T)
        }
        
        write(x = paste0("目录原始指向为:",filename," --> ",rawfile),
              file = linkinfo,
              append = T)
        
      }
      linkdirname <- paste0(filename,"-",format(Sys.time(), "%Y%m%d"),"-",Sys.info()["user"])
      
      if(file.exists(linkdirname)){
        linkdirname <- paste0(linkdirname,"-",RandomCode(n = 5))
      }
      
      createdir(filename = linkdirname,force = force)
      system(command = paste0("ln -s '",basename(linkdirname),"' '",filename,"'"))
      
      write(x = paste0("目录指向更改为:",filename," --> ",linkdirname),
            file = linkinfo,
            append = T)
      return(filename)
    }
    
    if (!file.exists(filename)) {
      print(paste0(filename, "文件夹进行创建"))
      
      if(!file.exists(dirname(filename))){
        
        createdir(filename = dirname(filename),
                  linkdir = F,
                  force = F)
        
      }
      dir.create(filename,
                 recursive = F,
                 mode = "0777")
      Sys.chmod(paths = filename,mode = "0777",use_umask = F)
      
    } else {
      # print(paste0(filename, "文件夹已存在"))
      
      if (force) {
        print(paste0(filename, "文件夹强制删除重新创建"))
        unlink(filename,
               recursive = T,
               force = T)
        dir.create(filename,
                   recursive = T,
                   mode = "0777")
        Sys.chmod(paths = filename,mode = "0777",use_umask = F)
      }
    }
    
  })
  
  return(filename)
}

#' 创建文件夹后进入文件夹
#'
#' @param filename 向量，文件夹名
#' @param force 逻辑，文件夹存在时是否强制删除重新创建，默认F
#' 
#' @return 向量,原始路径
#' 
#' @export
setwddir <- function(filename = "./",
                     force = F) {
  tryCatch({
    
    createdir(filename = filename,
              force = force)
    wd <- getwd()
    # print(paste0("移动到",filename,"文件夹"))
    setwd(filename)
    return(wd)
    
  })
  
}


#' 创建文件夹后进入文件夹运行函数
#'
#' @param path 向量，运行路径名
#' @param moudle 运行函数名
#' @param moudlename 运行函数名,字符串形式
#' @param force 逻辑，文件夹存在时是否强制删除重新创建，默认F
#' @param errorstop 逻辑，针对报错是否停止，默认为F
#' @param ... 传递进moudle
#' @export
runinpath <- function(path = "./",
                      moudle,
                      moudlename = "",
                      force = F,
                      errorstop = F,
                      simple = T,
                      ...) {
  print("--------------------")
  if(!simple){
    print("运行开始")
  }
  wd <- setwddir(filename = path, 
                 force = force)
  if(!simple){
    print(paste0("进入 ", path, " 位置"))
    starttime <- lubridate::now(tzone = "Asia/Shanghai")
    print(paste0("运行脚本:",moudlename,"   运行开始时间:",starttime))
  }
  
  data_try <- try({
    args <- moudle(...)
  },silent = F)
  
  if(!simple){
    endtime <- lubridate::now(tzone = "Asia/Shanghai")
    difftime <- format(difftime(endtime,starttime,units = "auto"),digits = 2)
    print(paste0("运行脚本:",moudlename,"   运行结束时间:",endtime,"   运行时间:",difftime))
    print(paste0("返回 ", wd, " 位置"))
  }
  
  setwd(wd)
  
  if ("try-error" %in% class(data_try) & errorstop) {
    stop(data_try)
  } else if ("try-error" %in% class(data_try)) {
    warning(data_try, immediate. = T)
    args <- NULL
  }
  
  if(!simple){
    print("运行结束")
  }
  return(args)
}
