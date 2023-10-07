#!/opt/conda/bin/Rscript

#' 随机森林分析
#'
#' @param filename 数据路径或数据
#' @param name 保存文件名
#' @param force 创建目录时是否清空原始目录
#' @param trainRatio 训练比例，1时全部为训练后做测试,改为0.7,后续可修改
#' @param importanceTop 选TOP的数量
#' @param rfPar 随机森林排序的方式，accu或gini
#' @param rep 折交叉检验重复数
#' @param cv.fold 折交叉检验次数
#' @param step 折交叉检验比例
#'
#' @export
randomForestAnalyst <- function(filename,
                                resultPath = "随机森林分析",
                                trainRatio = 0.7,
                                importanceTop = 30,
                                rfPar = "accu",
                                rep = 5,
                                cv.fold = 10,
                                step = 1.5,
                                num = NULL,
                                maxnum = 10,
                                imagetype = c("jpg", "pdf"),
                                ...){
  suppressMessages(library("randomForest"))
  suppressMessages(library("pROC"))
  suppressMessages(library("tidyverse"))
  suppressMessages(library("ggpubr"))
  
  wd <- getwd()
  
  data <- readdata(filename = filename, row.names = 1)
  
  # 创建结果目录
  setwddir(filename = resultPath,force = T)
  createdir("1.Importance")
  createdir("2.rfcv")
  createdir("3.SingleRoc")
  createdir("4.LogisticRoc")
  
  temp <- data[2:dim(data)[1],, drop = F]
  temp[,] <- apply(temp, 2, as.numeric)
  
  groupInfo <- data.frame(Sample = colnames(data),
                          Group = unlist(data[1,]))
  metaInfo <- data.frame("Name" = rownames(temp))
  metaInfo["name_new"] <- paste0("X",c(1:dim(metaInfo)[1]))
  
  data_all <- data.frame(as.matrix(t(temp)),check.names = F,stringsAsFactors = F)#数据由长变宽
  colnames(data_all) <- whto(metaInfo,colnames(data_all))
  data_all["sample"] <- rownames(data_all)
  data_all <- merge(groupInfo,data_all,by.x="Sample",by.y="sample",all.y=T)
  row.names(data_all) <- data_all[,"Sample"]
  data_all <- data_all[-1] #去除样本名称，留下分组信息
  
  #2.数据集拆分，默认7：3
  set.seed(1000)
  data_all$Group <- factor(data_all$Group)
  if(trainRatio == 1){
    data_train <- data_all
    data_test <- data_all
  }else{
    train <- sample(nrow(data_all), nrow(data_all)*trainRatio)
    data_train <- data_all[train, ]
    data_test <- data_all[-train, ]
  }
  
  #3.rf:代谢物重要性-----------------------------
  train_forest <- randomForest(Group~., data = data_train, importance = TRUE,ntree=(dim(data_train)[2]*5))
  train_importance <- data.frame(importance(train_forest), check.names = FALSE)
  train_importance <- train_importance[,c("MeanDecreaseAccuracy","MeanDecreaseGini")]
  train_importance["name_new"] <- rownames(train_importance)
  train_importance <- merge(train_importance,metaInfo,by.x="name_new",by.y="name_new",all.x=T)
  train_importance <- train_importance[order(train_importance$MeanDecreaseAccuracy, decreasing = TRUE), ]
  train_importance["Range"] <- c(1:(dim(train_importance)[1]))
  train_importance <- train_importance[c("Range","Name","MeanDecreaseAccuracy","MeanDecreaseGini")]
  savexlsx1(data = train_importance,
            filename = "./1.Importance/importance.xlsx",
            sheet = "importance")
  
  #importance绘图：
  #绘图时物质数量控制
  if(dim(train_importance)[1]>30 & (!is.na(importanceTop))){
    topN <- importanceTop
  }else{
    topN <- dim(train_importance)[1]
  }
  
  temp1 <- train_importance%>%arrange(desc(MeanDecreaseAccuracy))
  temp1 <- temp1[1:topN,]
  #避免ggplot根据字母表重排
  temp1$Name <- factor(temp1$Name,levels=rev(temp1$Name))
  p <- ggplot(data=temp1,mapping = aes(y=Name,x=MeanDecreaseAccuracy))+
    geom_point(shape=19)+
    theme_bw()+
    theme(axis.text=ggplot2::element_text(size=6),
          axis.title=ggplot2::element_text(size=10),
          legend.text=ggplot2::element_text(size=10),
          legend.title=ggplot2::element_text(size=10),
          aspect.ratio=16/9)+
    theme(legend.justification=c(1,0), legend.position=c(1,0))
  
  ggplotsave(plot = p,
             savepath = "1.Importance",
             mapname = "MeanDecreaseAccuracy", 
             imagetype = imagetype, 
             height = 4.5, width = 7,
             ...)
  
  temp2 <- train_importance%>%arrange(desc(MeanDecreaseGini))
  temp2 <- temp2[1:topN,]
  temp2$Name <- factor(temp2$Name,levels=rev(temp2$Name))
  p <- ggplot(data=temp2,mapping = aes(y = Name,x=MeanDecreaseGini))+
    geom_point(shape=19)+
    theme_bw()+
    theme(axis.text=ggplot2::element_text(size=6),
          axis.title=ggplot2::element_text(size=10),
          legend.text=ggplot2::element_text(size=10),
          legend.title=ggplot2::element_text(size=10),
          aspect.ratio=16/9)+
    theme(legend.justification=c(1,0), legend.position=c(1,0))
  
  ggplotsave(plot = p,
             savepath = "1.Importance",
             mapname = "MeanDecreaseGini", 
             imagetype = imagetype, 
             height = 4.5, width = 7,
             ...)
  
  
  #4.rfcv:折交叉验证---------------------------------
  train_cv <- replicate(rep,
                        rfcv(trainx = data_train[-ncol(data_train)],
                             trainy = data_train$Group,
                             cv.fold = cv.fold,
                             step = step),
                        simplify = FALSE)
  train_cv <- data.frame(sapply(train_cv, '[[', 'error.cv'))
  train_cv$datas <- rownames(train_cv)
  train_cv <- reshape2::melt(train_cv, id = 'datas')
  train_cv$datas <- as.numeric(as.character(train_cv$datas))
  cv_mean <- aggregate(train_cv$value, by = list(train_cv$datas), FUN = mean)
  colnames(cv_mean) <- c("number_of_datas","cross_validation_error")
  
  savexlsx1(data = cv_mean,
            filename = "./2.rfcv/cv_error.xlsx",
            sheet = "cv_error")
  
  p <- ggplot(cv_mean,aes(number_of_datas, cross_validation_error)) +
    geom_line() +
    theme_bw()+
    theme(aspect.ratio=9/16,
          panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(title = '',x = 'Number of datas', y = 'Cross-validation error')
  
  ggplotsave(plot = p,
             savepath = "2.rfcv",
             mapname = "cv_error", 
             imagetype = imagetype, 
             height = 4.5, width = 7,
             ...)
  
  #5.test set:单个代谢物ROC---------------------------
  # single
  imp <- max(which(cv_mean[2]==min(cv_mean)))
  # impMin <- min(which(cv_mean[2]==min(cv_mean)))+2  #多物质绘制多ROC曲线
  if(rfPar=="accu"){
    train_importance <- train_importance%>%arrange(desc(MeanDecreaseAccuracy))
  }else if(rfPar=="gini"){
    train_importance <- train_importance%>%arrange(desc(MeanDecreaseGini))
  }else{
    train_importance <- train_importance
  }
  # imp
  if(is.null(num)){
    imp <- cv_mean[imp,1]
  }else{
    imp <- num
  }
  
  if(is.na(imp)){
    imp <- maxnum
  }else if(imp >= maxnum){
    imp <- maxnum
  }else{
    imp <- imp
  }
  
  if(imp > dim(train_importance)[1]){
    imp <- dim(train_importance)[1]
  }
  
  metaname <- train_importance[1:imp,"Name"]
  data_test2 <- data_test
  colnames(data_test2)[-1] <- whto(metaInfo[,c(2,1)],colnames(data_test2)[-1])
  data_test2 <- as.data.frame(t(data_test2))
  data_test2 <- data_test2[c("Group",metaname),]
  
  # 3.SingleRoc
  auto_roc(data = data_test2,
           savepath = "3.SingleRoc",
           mapname = "Singleroc",
           imagetype = imagetype,
           ...)
  
  # 4.LogisticRoc
  auc <- NULL
  data2 <- NULL
  for (i in 1:imp) {
    aucdata2<- auto_logroc(data = data_test2[1:(i+1),],
                           savepath = "4.LogisticRoc",
                           mapname = paste0("LogisticRoc-",i),
                           imagetype = imagetype,
                           rocname = as.character(i),
                           saveroc = F,
                           ...)
    auc <- rbind(auc, aucdata2$auc)
    data2 <- rbind(data2, aucdata2$data)
  }
  savexlsx1(data = auc, filename = paste0("4.LogisticRoc/AUC.xlsx"), sheet = "AUC")
  savexlsx3(data = data2, filename = paste0("4.LogisticRoc/Logistic.xlsx"), sheet = "Logistic")
  
  data3 <- rbind(data_test2[1, ],data2)
  
  aucdata2 <- auto_roc(data = data3,
                       savepath = "4.LogisticRoc",
                       mapname = "OverlayROC",
                       imagetype = imagetype,
                       number = NULL,
                       saveroc = F,
                       height = 6,width = 6,
                       ...)
  
  p <- ggline(auc,
              x = "Feature", y = "AUC",
              xlab = "The Number of Feature",ylab = "AUC")+
    theme_bw()+
    theme(aspect.ratio=9/16,
          panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'))
  
  ggplotsave(plot = p,
             mapname = "AccumulationAUC",
             savepath = "4.LogisticRoc",
             imagetype = imagetype, 
             height = 4.5, width = 7,
             ...)
  
  # 生成报告
  mkreport(type = "随机森林")
  
  setwd(wd)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-f","--filename",default = "表达数据矩阵", nargs = "+",
                      help = "表达数据矩阵文件")
  parser$add_argument("-tr","--trainRatio", default = 0.7, type = "double", help="训练比例，1时全部为训练后做测试,改为0.7,后续可修改")
  parser$add_argument("-imp","--importanceTop", default = 30, type = "double", help="选TOP的数量")
  parser$add_argument("-rfp","--rfPar", default = "accu", help="随机森林排序的方式，accu或gini")
  parser$add_argument("-re","--rep", default = 5, type = "integer", help="折交叉检验重复数")
  parser$add_argument("-cv","--cv.fold", default = 10, type = "integer", help="折交叉检验次数")
  parser$add_argument("-st","--step", default = 1.5, type = "double", help=" 折交叉检验比例")
  parser$add_argument("-r","--resultPath",default = "随机森林分析", help = "结果输出路径")
  parser$add_argument("-z","--zip",default = F, help = "是否压缩",action='store_true')

  args <- parser$parse_args()
  
  zip <- args$zip
  args$zip <- NULL
  
  randomforestresults <- do.call(what = randomForestAnalyst, args = args)
  
  if(zip){
    zip::zip(zipfile = "分析结果.zip",files = args$resultPath)
    unlink(list.files(pattern = "[^z][^i][^p]$",include.dirs = T,all.files = T), recursive = T)
  }
}
