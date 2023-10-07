#! /opt/conda/bin/Rscript

#' @export
CV_plot <-function(
      inputfile = "./data.xlsx",
	  outpath = "./",
	  imagetype=c("png","pdf"),
	  width = 18,
      height = 15, 
      units = "cm",
      dpi=300
){
	library('ggplot2')
	library('dplyr')
	library(tidyverse)
	library(tidyr)
	library(lmbio)

	CumulativeFrequency <- list()
	freq_values <- c(0,0.05,0.1,0.15,0.2, 0.3,0.4, 0.5, 0.6,0.75, 1, 2)
	#qc rsd
	qc_data <- readxlsx(filename=inputfile,sheet="稳定性信息")%>%as.data.frame()
	qc_data_filtered <- qc_data[qc_data$`Ave(STD)` != 0, ]
	qc_rsd <- qc_data_filtered[["Rsd(STD)"]]/100
	qc_freq <- ecdf(qc_rsd)(freq_values)
	CumulativeFrequency[["QC"]]<- qc_freq

	rawdata<-readxlsx(filename=inputfile,sheet="定量信息")%>%as.data.frame()
	sample_group<-readxlsx(filename=inputfile,sheet="样本信息")%>%as.data.frame()#%>%column_to_rownames("样本分析名")
	#df_cv <- cbind(qc_rsd,sample_group)

	for ( group in unique(sample_group$"分组"[!grepl("STD-QC",sample_group$"分组")])){
	   samples <- sample_group[which(sample_group$"分组"==group),"样本分析名"]
	   group_data<- rawdata[,colnames(rawdata)%in% samples]
	   filtered <- subset(group_data, !rowSums(group_data == 0))
	   group_rsd <- apply(filtered,1,sd)/apply(filtered,1,mean)
	   group_freq<-ecdf(group_rsd)(freq_values)
	   CumulativeFrequency[[group]]<-group_freq
	}
	df_group_freq <- gather(data.frame(CumulativeFrequency), key = "Group", value = "Value")
	df_group_freq$freq <- rep(freq_values, length(unique(df_group_freq$Group)))
	#openxlsx::write.xlsx(df_group_freq,"df_group_freq.xlsx")


	pp = ggplot(df_group_freq, aes(x=as.numeric(freq), y = Value,color = Group))+
		geom_line() +
		#geom_smooth(method = "loess", se = FALSE) +  # 绘制曲滑的曲线
		theme_bw() + 
		theme(panel.grid.major = element_blank(), 
			  panel.grid.minor = element_blank(),
			  #panel.border = element_blank(),
			  axis.line=element_line(linewidth=0.6,colour="black"),
			  axis.text.x=element_text(colour="black",size = 8), 
			  axis.text.y=element_text(colour="black",size = 8), 
			  axis.title.x=element_text(size = 10), 
			  axis.title.y=element_text(size = 10)
		)+
		geom_vline(xintercept = 0.2, linetype = "dashed",colour = "black")+
		geom_vline(xintercept = 0.3,linetype = "dashed", colour = "black")+
		geom_hline(yintercept = 0.8, linetype = "dashed", color = "black")+
		scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1), limits = c(0, 1)) +  # 定义y轴刻度
		scale_x_continuous(breaks = c(0, 0.2,0.3,0.5,0.75,1,2)) +
		ylab("Cumulative Frequency") +
		xlab("CV")
		#scale_color_manual(values = c("QC" = "blue",unique(sample_group$"分组"[!grepl("STD-QC",sample_group$"分组")])))# 将QC的颜色设为蓝色
    for(j in 1:length(imagetype))
     {ggsave(file =paste0(outpath,"/CV_plot.",imagetype[j]),plot=pp,width = width, height = height, units = units,dpi=dpi)}

}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  parser$add_argument("-if","--inputfile", default = "./data.xlsx", help = "代谢数据文件，包括定量信息、样本信息和定量信息等sheet")
  parser$add_argument("-o","--outpath", default = "./", help = "输出路径")
  parser$add_argument("-w","--width",default = 18, help = "图片保存宽度")
  parser$add_argument("-he","--height",default = 15, help = "图片保存高度")
  parser$add_argument("-dpi","--dpi",default = 300, help = "分辨率")
  parser$add_argument("-u","--units",default = "cm", help = "基本单位")
  parser$add_argument("-it","--imagetype",default =c("png", "pdf"), help = "图片保存格式")
  args <- parser$parse_args()
  result <- do.call(what = CV_plot,args = args)
}