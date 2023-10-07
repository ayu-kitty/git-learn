#!/opt/conda/bin/Rscript
#' @export
frequencydraw <- function(file = "shift_all_cor_rsd.csv",
                          savepath = "./",
                          filename = "Cumulative Frequency",
                          imagetype = c("jpg","pdf"),
                          point_shape = 17,
                          line_colour = "red",
                          line_size = 0.5,
                          width = 10,
                          height = 8, 
                          family = "sans",
                          vline_yintercept = 30,
                          vline_lty = 1,
                          vline_colour = "red",
                          ...){
  
  suppressWarnings(library(ggplot2))
  
  data <- readdata(file)
  # data <- data[!colnames(data) %in% "adduct"]
  # data <- subset(data, data$name!="")
  QC_curve <- dplyr::select(data,contains("qc_rsd",ignore.case = T))
  pdata <- table(QC_curve$qc_rsd)
  pframe <- data.frame(as.numeric(names(pdata)),as.numeric(pdata))
  colnames(pframe) <- c("qc_rsd","num")
  rankData <- pframe[order(pframe[,1],decreasing=F),]
  mulData <- data.frame()
  sum <- 0
  for(i in 1:nrow(rankData)){
    sum <- sum+rankData[i,2]
    row <- cbind(rankData$qc_rsd[i],sum)
    mulData <- rbind(mulData,row)
  }
  colnames(mulData) <- c("qc_rsd","num")
  curve_Data  <- data.frame(RSD =mulData$qc_rsd,PER =mulData$num/sum*100)
  pc <- ggplot(curve_Data, aes(x = RSD, y = PER))+
    geom_point(shape=point_shape) +
    geom_line(colour = line_colour,size = line_size) +
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line=element_line(size=0.4,colour="black"),
          axis.text.x=element_text(colour="black",size = 6), 
          axis.text.y=element_text(colour="black",size = 6), 
          axis.title.x=element_text(size = 8), 
          axis.title.y=element_text(size = 5))+
    geom_vline(xintercept = vline_yintercept,lty =vline_lty, colour = vline_colour)+
    ylab("cumulative frequency") +
    xlab("RSD (%)") 
  
  ggplotsave(plot = pc,
             savepath = savepath,
             mapname = filename,
             imagetype = imagetype,
             width = width, 
             height = height,
             family = family,
             ...)
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  library("argparse")
  parser <- ArgumentParser()
  parser$add_argument("-f","--file",default = "shift_all_cor_rsd.csv", help = "矫正后数据矩阵具有rsd信息")
  parser$add_argument("-sp","--savepath",default = "./", help = "保存数据路径")
  parser$add_argument("-fn","--filename",default = "Cumulative Frequency", help = "保存图片名称")
  parser$add_argument("-it","--imagetype",default =c("jpg", "pdf"), help = "图片保存格式",nargs= "+")
  parser$add_argument("-ps","--point_shape",default = 17, help = "点的形状")
  parser$add_argument("-lc","--line_colour",default = "red", help = "线条颜色")
  parser$add_argument("-ls","--line_size",default = 0.5, help = "线条粗细")
  parser$add_argument("-w","--width",default = 12, help = "图片保存宽度")
  parser$add_argument("-he","--height",default = 8, help = "图片保存高度")
  parser$add_argument("-dpi","--dpi",default = 300, help = "分辨率")
  parser$add_argument("-u","--units",default = "cm", help = "基本单位")
  parser$add_argument("-fa","--family",default = "sans", help = "字体样式")
  
  args <- parser$parse_args()
  result <- do.call(what = frequencydraw,args = args) 
} 
