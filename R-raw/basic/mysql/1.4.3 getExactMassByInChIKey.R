#!/opt/conda/bin/Rscript

#' 根据inchikey查询精确分子量
#' 
#' @param InChIKey 单个inchikey 或者 向量Vector<inchikey>
#'
#' @return dataframe
#'                      inchikey  exactMass
#' 1 RAKAWZJMWJKWRH-QRMJXLNNSA-N 684.565253
#' 2 VOGBKCAANIAXCI-UHFFFAOYSA-N 626.491025
#' ...
#' 
#' @export
getExactMassByInChIKey <- function(InChIKey){
  con <- linkmysql(dbname = "cosa")
  tmp <- "SELECT b.inchikey,a.exactMass FROM compound_property AS a INNER JOIN compound_structure AS b ON a.cid = b.cid WHERE b.inchikey"
  sql <- paste0(tmp, " IN (\"", paste(InChIKey,collapse = "\",\""),"\")" )
  rs <- dbSendStatement(con,sql)
  data <- dbFetch(rs,n = -1)
  dbClearResult(rs)
  dbDisconnect(con)
  return(data)
}
