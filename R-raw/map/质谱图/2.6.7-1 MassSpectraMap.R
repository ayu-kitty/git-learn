#!/opt/conda/bin/Rscript

# MChromatograms
# OnDiskMSnExp

MSplot <- function(x,
                   y,
                   group = NULL,
                   breaks = seq(0, 45, 5),
                   showrt = F,
                   distance = 20,
                   size = 0.5,
                   color = "#FF6600",
                   mapname = "MassSpectra",
                   width = 8,
                   height = 4,
                   savepath = "./",
                   ...) {
  suppressMessages(library("ggplot2"))
  suppressMessages(library("ggrepel"))
  # x <- x/60
  if (is.null(group)) {
    data <- data.frame(`Retentime time(min)` = x,
                       Intensity = y,
                       check.names = F)
  } else {
    data <- data.frame(`Retentime time(min)` = x,
                       Intensity = y,
                       Sample = group,
                       check.names = F)
  }
  
  data <- data[order(data$`Retentime time(min)`), ]
  data1 <- data.frame()
  
  mintime <- min(data$`Retentime time(min)`)
  maxtime <- max(data$`Retentime time(min)`)
  
  t <- mintime
  
  while (t < maxtime) {
    data2 <- data[data$`Retentime time(min)` >= t, ]
    t <- t + (maxtime - mintime) / 1000
    data2 <- data2[data2$`Retentime time(min)` < t, ]
    if (is.null(group)) {
      data2 <- data2[which.max(data2$Intensity), ]
      data1 <- rbind(data1, data2)
    } else {
      for (group1 in unique(group)) {
        data3 <- data2[data2$Sample == group1, ]
        data3 <- data3[which.max(data3$Intensity), ]
        data1 <- rbind(data1, data3)
      }
    }
  }
  
  data <- data1
  data[, "tag"] <- ""
  mintime <- min(data$`Retentime time(min)`)
  maxtime <- max(data$`Retentime time(min)`)
  maxIntensity <- max(data$Intensity)
  
  t <- mintime
  
  while (t < maxtime) {
    data2 <- data[data$`Retentime time(min)` >= t, ]
    t <- t + (maxtime - mintime) / distance
    data2 <- data2[data2$`Retentime time(min)` < t, ]
    data2 <- data2[which.max(data2$Intensity), ]
    
    if (1 %in% which(data$`Retentime time(min)` == data2$`Retentime time(min)`)) {
      
    } else if (dim(data)[1] %in% which(data$`Retentime time(min)` == data2$`Retentime time(min)`)) {
      
    } else {
      if ((data[which(data$`Retentime time(min)` == data2$`Retentime time(min)`)[1] - 1, "Intensity"] < data2$Intensity) &
          (data[which(data$`Retentime time(min)` == data2$`Retentime time(min)`)[1] + 1, "Intensity"] <= data2$Intensity)) {
        data[data$`Retentime time(min)` == data2$`Retentime time(min)`, "tag"] <- ifelse(data2$Intensity >= (maxIntensity / 20), 
                                                                                         round(data2$`Retentime time(min)`, 4), 
                                                                                         "")
      }
    }
  }
  
  if (is.null(group)) {
    p <- ggplot(data, aes(x = `Retentime time(min)`, y = Intensity)) +
      geom_line(color = color, size = size)
  } else {
    p <- ggplot(data, aes(x = `Retentime time(min)`, y = Intensity)) +
      geom_line(mapping = aes(color = Sample), size = size)
  }
  
  p <- p +xlab("Retention time(min)") +
    ylab("Intensity") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))

  # LC项目,最大maxtime 在15~20内,只保留15min的
  if (maxtime > 15 & maxtime <20){
    p <- p + scale_x_continuous(breaks = breaks,limits = c(0,15))
  }else{
    p <- p + scale_x_continuous(breaks = breaks,limits = c(0,maxtime))
  }
  
  if (showrt) {
    p <- p + geom_text_repel(aes(label = tag), size = 2.5, direction = "y")
  }
  
  ggplotsave(plot = p,
             mapname = mapname,
             height = height,
             width = width,
             savepath = savepath,
             ...)
  
}

#' @export
mzplot <- function(data,
                   ...) {
  UseMethod("mzplot")
}

#' @export
mzplot.default <- function(data,...){
  print("~无相应类型数据类型进行质谱图绘制")
}

#' @export
mzplot.MSnExp <- function(data,
                          type = "BPC",
                          gap = 0.05,
                          ...) {
  
  if (type == "BPC") {
    msdata <- MSnbase::chromatogram(data, aggregationFun = "max")
  } else if (type == "TIC") {
    msdata <- MSnbase::chromatogram(data, aggregationFun = "sum")
  }
  
  mapdata <- NULL
  
  for ( i in 1:length(msdata)) {
    mapdata1 <- data.frame(rtime = rtime(msdata[[i]])/60,
                           intensity = intensity(msdata[[i]]),
                           group = colnames(msdata)[i],
                           check.names = F,
                           stringsAsFactors = F)
    mapdata <- rbind(mapdata,mapdata1)
  }
  
  n <- 0
  idata <- max(mapdata$intensity,na.rm = T)
  for (group1 in unique(mapdata$group)) {
    mapdata[mapdata$group == group1,"intensity"] <- mapdata[mapdata$group == group1,"intensity"]+idata*(gap*n)
    n <- n+1
  }
  
  if(length(msdata) > 1){
    MSplot(x = mapdata$rtime,y = mapdata$intensity,group = mapdata$group,...)
  }else{
    MSplot(x = mapdata$rtime,y = mapdata$intensity,...)
  }
  
}

#' @export
mzplot.MChromatograms <- function(data,
                                  gap = 0.05,
                                  type = NULL,
                                  ...) {
  mapdata <- NULL
  
  # data2 <<- data
  
  for ( i in 1:dim(data)[2]) {
    for ( j in 1:dim(data)[1]) {
      if(length(rtime(data[j,i]))==0){
        next
      }
      mapdata1 <- data.frame(rtime = rtime(data[j,i]),
                             intensity = intensity(data[j,i]),
                             group = colnames(data)[i],
                             check.names = F,
                             stringsAsFactors = F)
      mapdata <- rbind(mapdata,mapdata1)
    }
  }
  
  n <- 0
  idata <- max(mapdata$intensity,na.rm = T)
  for (group1 in unique(mapdata$group)) {
    mapdata[mapdata$group == group1,"intensity"] <- mapdata[mapdata$group == group1,"intensity"]+idata*(gap*n)
    n <- n+1
  }
  
  if(dim(data)[2] > 1){
    MSplot(x = mapdata$rtime,y = mapdata$intensity,group = mapdata$group,...)
  }else{
    MSplot(x = mapdata$rtime,y = mapdata$intensity,...)
  }
}

#' @export
mzplot.OnDiskMSnExp <- mzplot.MSnExp

#' @export
map_mass_massmap <- function(filename,
                             samplename = NULL,
                             ...){
  suppressMessages(library("MSnbase"))
  
  data <- readdata(filename = filename,samplename = samplename)
  
  result <- mzplot(data,...)
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "mzML格式文件",
                      required = T)
  parser$add_argument("-sa","--samplename",default = NULL,nargs="+",
                      help = "mzML文件对应样本名")
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("png","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf","html"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-g","--gap",default = 0.05,type= "double",help = "多样本间隔")
  parser$add_argument("-sr","--showrt", default = F,action = "store_true", help="是否展示保留时间")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  result <- do.call(what = map_mass_massmap, args = args)
}
