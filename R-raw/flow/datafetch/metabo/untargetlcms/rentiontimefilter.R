#' @export
rentiontimefilter <- function(data,
                              rttimemin = 0.5,
                              rttimemax= 15.5){
  
  data1 <- data
  dataraw<-data1$data$rawdata$data
  filternum<-which(dataraw$`Retention time (min)`<=rttimemax&dataraw$`Retention time (min)`>=rttimemin)
  rtfilter<-dataraw[filternum,]
  data$data$rawdata$data<-rtfilter
  print("数据筛选保留时间完毕")
  return(data)
  
}
