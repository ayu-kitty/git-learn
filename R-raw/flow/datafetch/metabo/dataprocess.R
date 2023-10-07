#' autodataprocess
#'
#' 数据获取并进行预处理
#'
#' @param data obj
#' @param fetch 是否数据获取
#' @param name 数据获取xlsx名称，不填写将从data中获取路径
#' @param sheet xlsx的sheet名
#' @param order 样本名排序
#' @param del 逻辑，是否删除
#' @param mean 逻辑，是否均值
#' @param missvalue 0值筛选标准
#' @param group 逻辑，是否根据组别进行0值筛选
#' @param rsd rsd筛选标准
#' @param zeroprocess 0值处理方法
#' @param ...
#'
#' @export
autodataprocess <- function(data,
                            fetch = T,
                            name = NA,
                            sheet = 1,
                            order = T, del = T, mean = T,
                            missvalue = data$info$handle$missvalue, group = T,
                            rsd = data$info$handle$rsd,
                            zeroprocess = data$info$handle$zeroprocess,
                            saminfo = "批次信息.xlsx",
                            method = "statTarget",
                            minqc = 0,
                            minsample = 0,
                            Frule = 0.8, 
                            MLmethod = "QCRFSC", 
                            QCspan = 0,
                            imputeM = "KNN",
                            adjustsavepath = "批次校正",
                            ...) {
  print("********数据矩阵处理开始********")
  
  if (fetch) {
    if (!is.na(name)) {
      data$data$rawdata$data <- readdata(name = name, sheet = sheet)
    } else {
      data <- datafetch(data, ...)
    }
  }
  
  if(file.exists(saminfo)){
    data$data$rawdata$data <- Batch_correct(input = data$data$rawdata$data,
                                            saminfo = saminfo,
                                            method = method,
                                            minqc = minqc,
                                            minsample = minsample,
                                            Frule = Frule, 
                                            MLmethod = MLmethod, 
                                            QCspan = QCspan,
                                            imputeM = imputeM,
                                            savepath = adjustsavepath)
  }
  
  # if(data[["info"]][["basic"]][["项目类型"]]=="非靶向代谢-LCMS"||data[["info"]][["basic"]][["项目类型"]]=="非靶向代谢-LCMS-EMBD"||data[["info"]][["basic"]][["项目类型"]]=="非靶向代谢-LCMS-PMDB")
  # {data<-rentiontimefilter(data,rttimemin=rttimemin,rttimemax=rttimemax)}
  
  data <- preprocess(data, order = order, del = del, mean = mean)
  # data <- annotatekegg(data)
  data <- automissvalue(data, missvalue = missvalue, group = group)
  data <- autorsdfilter(data, rsd = rsd)
  data <- zeroprocess(data, zeroprocess = zeroprocess)
  
  print("********数据矩阵处理结束********")
  return(data)
}
