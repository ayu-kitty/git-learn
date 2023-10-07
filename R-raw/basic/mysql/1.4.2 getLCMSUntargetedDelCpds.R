#!/opt/conda/bin/Rscript

#' 获取非靶LC-MS项目中要剔除物质的inchikey
#' 
#' @return dataframe
#'      cid                    inchikey
#' 1    47059 VEXZGXHMUGYJMC-UHFFFAOYSA-M
#' 2    47085 JLVVSXFLKOJNIY-UHFFFAOYSA-N
#'
#' @export
getLCMSUntargetedDelCpds <- function(){
  con <- linkmysql(dbname = "cosa")
  sql <- "SELECT t1.cid,t2.inchikey FROM (SELECT a.cid FROM compound_classification AS a WHERE a.kingdom = 'Inorganic compounds' UNION SELECT b.cid FROM compound_annotation AS b ) AS t1 INNER JOIN compound_structure AS t2 ON t1.cid = t2.cid"
  rs <- dbSendStatement(con,sql)
  data <- dbFetch(rs,n = -1)
  dbClearResult(rs)
  dbDisconnect(con)
  return(data)
}


#' 获取非靶LC-MS项目中要剔除的物质的hmdb编号
#' 
#' @export
getLCMSUntargetedDelCpdsforhmdbid <- function(table = "compound_identifier",
                                              id = "hmdbid"){
  delcpd <- getLCMSUntargetedDelCpds()
  
  delcpd <- getmysqldata(dbname = "cosa",
                         table = table,
                         wherename = "cid",wheredata = delcpd$cid)
  
  delcpdid <- delcpd[,id]
  delcpdid <- delcpdid[!is.na(delcpdid)]
  delcpdid <- delcpdid[!duplicated(delcpdid)]
  
  return(delcpdid)
}


#' 获取非靶LC-MS项目中要剔除的物质的inchikey编号
#' 
#' @export
getLCMSUntargetedDelCpdsforinchikey <- function(table = "compound_structure",
                                                id = "inchikey"){
  delcpdid <- getLCMSUntargetedDelCpdsforhmdbid(table = table,
                                                id = id)
  
  return(delcpdid)
}