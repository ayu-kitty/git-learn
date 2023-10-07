#!/opt/conda/bin/Rscript

#' 生成空代项目报告
#'
#' @param mode 报告模板
#' @param type 空代产品线类型,默认为空间代谢组(or 中药空间代谢组)
#' @export

makespacereport <- function(mode = "LM",
                            type = "空间代谢组",
							CLOUD=F
) {
  
  #7.1日修改报告名称为订单号-客户名称-联系人名称-项目类型项目报告
  report <- getreportname()
  
  setwd(report)
  
  mkreport(type = type,CLOUD=CLOUD)
  
  setwd("../")
  
}

