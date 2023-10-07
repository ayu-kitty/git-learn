#!/opt/conda/bin/Rscript

#' 获取项目报告名称
#'
#' @param filename 项目登记单路径
#' @param sheet 读取sheet
#'
#' @export
getreportname <- function(filename = "项目登记单.xlsx",
                          sheet = "项目登记单"){
  data <- readdata(filename = filename,
                   sheet = sheet)
  
  #7.1日修改报告名称为订单号-客户名称-联系人名称-项目类型项目报告
  if(!is.na(data[data[, 1] == "客户名称" & !is.na(data[, 1]), 2])){
    if(data[data[, 1] == "客户名称" & !is.na(data[, 1]), 2] == data[data[, 1] == "联系人名称" & !is.na(data[, 1]), 2]){
      
      report <- paste0(data[data[, 1] == "项目编号", 2],
                       "-", data[data[, 1] == "客户名称", 2],
                       "-", data[data[, 1] == "项目类别", 2], "(项目报告)")[1]
      
    }else{
      
      report <- paste0(data[data[, 1] == "项目编号", 2],
                       "-", data[data[, 1] == "客户名称", 2],
                       "-", data[data[, 1] == "联系人名称", 2],
                       "-", data[data[, 1] == "项目类别", 2], "(项目报告)")[1]
      
    }
  }else{
    report <- "项目报告"
  }
  
  return(report)
}
