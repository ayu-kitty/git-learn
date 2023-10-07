#!/opt/conda/bin/Rscript

#' 连接mysql数据库
#'
#' @param dbname 数据库名
#'
#' @export
linkmysql <- function(dbname = "meta"){
  suppressMessages(library("RMySQL"))
  con <- dbConnect(MySQL(),
                   user = 'lujw',
                   password = 'lumingbio',
                   dbname = dbname,
                   port = 3306,
                   host = '10.100.10.42')
  dbSendQuery(con, "SET NAMES utf8mb4")
  return(con)
}
