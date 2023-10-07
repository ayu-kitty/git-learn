
#' 获取vip
#'
#' @param data obj
#'
#' @export
getvip <- function(data) {
  
  if(!is.null(data$statistics)){
    mode <- data$statistics@typeC
    if(mode != "PCA"){
      vipdata <- data$statistics@vipVn
      vipdata <- data.frame(VIP = vipdata,
                            row.names = names(vipdata),
                            check.names = F,
                            stringsAsFactors = F)
      vipdata <- list(vipdata = vipdata,
                      mode = mode)
    }else{
      vipdata <- NULL
    }
  }else{
    vipdata <- NULL
  }
  
  return(vipdata)
}