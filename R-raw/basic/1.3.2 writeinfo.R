#!/opt/conda/bin/Rscript

#' 输出运行信息
#'
#' @param command 运行命令
#' @param filename 保存文件名
#' @param version 是否输出版本号
#' @param user 是否输出运行用户名
#' @param starttime 是否输出开始时间
#' @param endtime 是否输出结束时间
#'
#' @export
writeinfo <- function(command = commandArgs(),
                      filename = "runinfo.txt",
                      version = T,
                      user = T,
                      image = T,
                      starttime = T,
                      endtime = F){
  command <- paste(command,collapse = " ")
  command <- paste0("---------\n运行命令:",command)
  if(version){
    pkgname <- Biobase::package.version("lmbio")
    command <- paste0(command,"\n运行版本:lmbio-",pkgname)
  }
  
  if(user){
    command <- paste0(command,"\n运行用户:",Sys.info()["user"])
  }
  
  if(image){
    command <- paste0(command,"\n运行镜像:",system("source /etc/profile;echo $IMAGE_NAME",intern = T))
  }
  
  if(starttime){
    time <- Sys.time()
    command <- paste0(command,"\n开始时间:",time)
  }
  
  if(endtime){
    time <- Sys.time()
    command <- paste0("结束时间:",time)
  }
  
  write(x = command,file = filename,append = T)
  # return(command)
}
