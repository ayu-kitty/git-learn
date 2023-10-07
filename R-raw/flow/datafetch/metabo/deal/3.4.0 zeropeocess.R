#' @export
zeroprocess <- function(data,
                        zeroprocess = data$info$handle$zeroprocess) {
  data3 <- data
  data3$info$handle$zeroprocess <- zeroprocess


  if (!is.na(data3$info$handle$zeroprocess)) {
    print("0值处理开始")
	prezerodata<-data3$data$data$data
    data1 <- data3$data$data$data
     
    if (data3$info$handle$zeroprocess == "min") {
      print("半值处理模式")
      data1[data1 == 0] <- 0.5 * min(data1[data1 != 0])
    } else if (data3$info$handle$zeroprocess == "knn") {
      data1 <- t(data1)
      data1[data1 == 0] <- NA
      khan.imputed <- impute::impute.knn(as.matrix(data1))
      data1 <- t(khan.imputed[["data"]])
    } else {
      print("模式选择错误，仅有min，knn")
    }

    data3$data$zeroprocess$data <- data3$data$data$data
    data3$data$zeroprocess$dealdata <- data1
	data3$data$zeroprocess$prezerodata<-prezerodata
    data3$data$data$data <- data1
    print("0值处理完毕")
  } else {
    warning("未进行0值处理，请确认", immediate. = T)
  }

  return(data3)
}
