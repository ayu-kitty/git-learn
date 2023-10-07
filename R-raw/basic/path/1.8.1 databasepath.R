
#' 根据系统选择系统数据位置
#'
#' @param path 向量，路径
#'
#' @return 返回路径
#' @export
ospath <- function(path = NULL,exist = T) {
  
  if(!file.exists(path) & exist){
    stop(paste0(path,"文件不存在"))
  }
  
  return(path)
}


#' 数据库位置
#'
#' @param database 向量，数据库
#' @param path 向量，路径
#'
#' @return 返回数据库中路径
#' @export
databasepath <- function(database = "database/",
                         path = "",
                         exist = T) {
  if(path == ""){
    path1 <- paste0(lmbio::path$database,"/",database)
  }else{
    path1 <- paste0(lmbio::path$database,"/",database,"/",path)
  }
  
  path1 <- ospath(path1,exist = exist)
  
  return(path1)
}

#' 库位置
#'
#' @param path 向量，路径
#'
#' @return 返回数据库中路径
#' @export
packagepath <- function(path = "") {
  path1 <- paste0(lmbio::path$package,"/",path)
  path1 <- ospath(path1)
  
  return(path1)
}
