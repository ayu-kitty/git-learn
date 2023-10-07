#!/opt/conda/bin/Rscript

#' @export
Samplecorrplot2<-function(input,
                          type=NULL,
                          saminfo,
                          savepath="./"){
  
  suppressMessages(library('dplyr'))
  suppressMessages(library('openxlsx'))
  suppressMessages(library('ggplot2'))
  suppressMessages(library('corrplot'))
  suppressMessages(library('reshape'))
  suppressMessages(library('readxl'))
  #legend颜色
  if (is.null(type)){
    col11=c("#00468BFF","#0099B4FF","cornsilk3","#f19b78","#E64B35FF")
  }else{
    col11 <- read.xlsx("color_sample.xlsx",sheet=1)
    col11 <- col11[,1]
  }
  plot_path <- paste0(savepath,"./Samplecorr/")
  if (dir.exists(plot_path)){
    list.files("Samplecorr")
  }else{
    dir.create(plot_path)
  }
  #表达矩阵dataframe
  data<-input
  #读取分组信息,两列:第一列samle,第二列class
  sample_infor <- saminfo
  sample_infor<-sample_infor[!duplicated(sample_infor[,1]),]
  group <- unique(sample_infor[,2])#组名
  sample <- sample_infor[,1]#样本名
  targetdata1 <- data[,sample]
  #给每组赋予颜色
  if(!"color" %in% colnames(sample_infor)){
    groupColors <- setNames(SelectColors("visium9",n = length(group)), group)
    sample_infor<- sample_infor %>%mutate(color = groupColors[class])
  }
  ##参数设置
  sample_number <- c(2:200)
  number0 <- seq(48,0.5,length.out = 199)
  number0[6:9] <- number0[6:9]-8
  number0[4:5] <- number0[4:5]-15
  number0[2:3] <- number0[2:3]-20
  number0[1] <- number0[1]-45
  tl.cex <- number0[which(length(sample)==sample_number)]/length(sample)
  clcex=2
  #空值最小值填充
  targetdata1[is.na(targetdata1)] <- min(targetdata1,na.rm = T)
  savedata<-as.data.frame(cbind(rownames(data),targetdata1))
  names(savedata)[1]<-"name"
  write.xlsx(savedata,paste0(plot_path,"/Samplecorrdata.xlsx"))
  corr_data <- as.matrix(cor(targetdata1,method = "pearson"))
  res1 <- cor.mtest(corr_data,conf.level = .95)
  clcex=2
  if (length(sample)>100 & length(sample)<=200) {
    plotfile(
      savepath = plot_path,
      mapname = "Samplecorrplot",
      imagetype = c("jpg", "pdf"),
      height = 30,
      width = 30,
      dpi = 300,
      family = "sans",
      units = "in"
    )
    pp=corrplot::corrplot(corr_data,type="upper",tl.pos="lt",order="original",
                          insig = "label_sig",
                          sig.level = c(.001,.01,.05),
                          col=col11,
                          method = "circle",
                          tl.cex = tl.cex/(5/4),cl.cex = clcex,
                          tl.col = sample_infor$color)
    pp=corrplot::corrplot(corr_data,add=TRUE,type="lower", method="number",
                          order="original",diag=F,
                          col=col11,
                          number.cex = tl.cex/(5/4),cl.cex = clcex,
                          cl.pos = "n",tl.pos = "n")
    plotsave()
  }else if (length(sample)>200) {
    plotfile(
      savepath = plot_path,
      mapname = "Samplecorrplot",
      imagetype = c("jpg", "pdf"),
      height = 30,
      width = 30,
      dpi = 300,
      family = "sans",
      units = "in"
    )
    pp=corrplot::corrplot(corr_data,type="upper",tl.pos="lt",order="original",
                          insig = "label_sig",
                          sig.level = c(.001,.01,.05),
                          col=col11,
                          method = "circle",
                          tl.cex = 0.01,cl.cex = clcex,
                          tl.col = sample_infor$color)
    pp=corrplot::corrplot(corr_data,add=TRUE,type="lower", method="number",
                          order="original",diag=F,
                          col=col11,
                          number.cex = 0.01,cl.cex = clcex,
                          cl.pos = "n",tl.pos = "n")
    plotsave()
  }else{
    plotfile(
      savepath = plot_path,
      mapname = "Samplecorrplot",
      imagetype = c("jpg", "pdf"),
      height = 30,
      width = 30,
      dpi = 300,
      family = "sans",
      units = "in"
    )
    pp=corrplot::corrplot(corr_data,type="upper",tl.pos="lt",order="original",
                          p.mat = res1$p,
                          insig = "label_sig",
                          sig.level = c(.001,.01,.05),
                          pch.col = "white",
                          pch.cex = tl.cex/(6/4),
                          col=col11,
                          method = "circle",
                          tl.cex = tl.cex/(5/4),cl.cex =clcex,
                          tl.col = sample_infor$color)
    pp=corrplot::corrplot(corr_data,add=TRUE,type="lower", method="number",
                          order="original",diag=F,
                          col=col11,
                          number.cex = tl.cex/(6/4),cl.cex = clcex,
                          cl.pos = "n",tl.pos = "n")
    plotsave()
  }
  wb<-createWorkbook("./Samplecorr/Samplecorrplot.xlsx")
  addWorksheet(wb, "相关性系数")
  writeData(wb, sheet = "相关性系数", corr_data, rowNames = T, borders = "columns")
  addWorksheet(wb, "相关性检验")
  writeData(wb, sheet = "相关性检验", res1$p, rowNames = T, borders = "columns")
  saveWorkbook(wb, paste0(savepath,"./Samplecorr/Samplecorrplot.xlsx"), overwrite = T)
}
