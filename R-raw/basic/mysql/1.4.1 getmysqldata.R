#!/opt/conda/bin/Rscript

#' 获取mysql数据库的数据
#' 
#' @param dbname 数据库名
#' @param table 表格名
#' @param tablelist 表格列名
#' @param wherename 筛选列名
#' @param wheredata 筛选数据名
#' @param n 最多数量
#' @param limit 限制
#' @param ... 
#'
#' @export
getmysqldata <- function(dbname = "meta",
                         table = "adduct",
                         tablelist = "*",
                         wherename = NULL,
                         wheredata = NULL,
                         wherenotin = F,
                         n = -1,
                         limit = NULL,
                         ...){
  options(scipen = 10)
  con <- linkmysql(dbname = dbname)
  
  try({
    if(is.null(wherename)){
      Statement <- paste0("SELECT ",tablelist," FROM ",table)
    }else{
      if(wherenotin){
        Statement <- paste0("SELECT ",tablelist," FROM ",table," WHERE ",wherename," NOT IN (\"",paste(wheredata,collapse = "\",\""),"\")")
      }else{
        Statement <- paste0("SELECT ",tablelist," FROM ",table," WHERE ",wherename," IN (\"",paste(wheredata,collapse = "\",\""),"\")")
      }
    }
    
    if(!is.null(limit)){
      Statement <- paste0(Statement," LIMIT ",limit[1],",",limit[2])
    }
    
    # print(Statement)
    rs <- dbSendStatement(con,
                          Statement)
    
    data <- dbFetch(rs,n = n)
    data[data == ""] <- NA
    dbClearResult(rs)
  },silent = F)
  
  dbDisconnect(con)
  
  return(data)
}

#' getmysqlrow
#' 
#' 获取表格行数
#' 
#' @param dbname 数据库名
#' @param table 表格名
#'
#' @export
getmysqlrow <- function(dbname = "meta",
                        table = "adduct"){
  con <- linkmysql(dbname = dbname)
  
  try({
    rs <- dbSendStatement(
      con,
      paste0("SELECT COUNT(*) FROM ",table)
    )
    
    data <- dbFetch(rs)
    numrow <- data[1,1]
  },silent = F)
  
  dbDisconnect(con)
  return(numrow)
}

#' 模糊获取mysql数据库的数据
#' 
#' @param dbname 数据库名
#' @param table 表格名
#' @param tablelist 表格列名
#' @param wherename 筛选列名
#' @param wheredata 筛选数据名
#' @param n 最多数量
#' @param limit 限制
#' @param ... 
#'
#' @export
getmysqldata_vague <- function(dbname = "meta",
                               table = "adduct",
                               tablelist = "*",
                               wherename = NULL,
                               wheredata = NULL,
                               n = -1,
                               limit = NULL,
                               ...){
  options(scipen = 10)
  con <- linkmysql(dbname = dbname)
  
  try({
    if(is.null(wherename)){
      Statement <- paste0("SELECT ",tablelist," FROM ",table)
    }else{
        Statement <- paste0("SELECT ",tablelist," FROM ",table," WHERE REGEXP_LIKE(",wherename,",'(",paste0(wheredata,collapse = "|"),")')")
    }
    
    if(!is.null(limit)){
      Statement <- paste0(Statement," LIMIT ",limit[1],",",limit[2])
    }
    
    # print(Statement)
    rs <- dbSendStatement(
      con,
      Statement
    )
    
    data <- dbFetch(rs,n = n)
    data[data == ""] <- NA
    dbClearResult(rs)
  },silent = F)
  
  dbDisconnect(con)
  
  return(data)
}
