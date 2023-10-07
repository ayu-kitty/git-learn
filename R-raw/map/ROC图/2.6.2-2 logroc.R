#!/opt/conda/bin/Rscript

#' 逻辑回归后roc可视化
#'
#' @param data 数据
#' @param name 保存文件名
#' @param type 保存图片格式
#' @param saveroc 保存roc数据
#' @param dealname 是否处理行名
#' @param rocname 图像名称
#' @param ... 见[auto_plotroc()]
#'
#' @export
auto_logroc <- function(data,
                        mapname = "LogsticROC",
                        saveroc = T,
                        rocname = F,
                        savepath = "./",
                        ...) {
  
  data1 <- data
  
  group <- table(as.vector(t(data1[1, ])))
  
  if (any(group < 3)) {
    savetxt(data = "有比较组分析未提供roc分析，可能由于分组样本小于3个或比较组大于两组",
            filename = paste0(savepath,"/说明.txt"))
    warning("有分组样本小于3个", immediate. = T)
    return("有分组样本小于3个")
  } else if (length(group) > 2) {
    savetxt(data = "有比较组分析未提供roc分析，可能由于分组样本小于3个或比较组大于两组",
            filename = paste0(savepath,"/说明.txt"))
    warning("比较组大于两组", immediate. = T)
    return("比较组大于两组")
  }
  
  data1 <- as.data.frame(t(data1))
  data1[, 2:dim(data1)[2]] <- apply(data1[, 2:dim(data1)[2], drop = F], 2, as.numeric)
  
  if(rocname == F){
    metaname <- paste(colnames(data1)[-1], collapse = " + ")
  }else{
    metaname <- rocname
    # print(metaname)
  }
  
  colnames(data1)[-1] <- paste0("X", 1:(ncol(data1) - 1))
  formula <- paste(colnames(data1)[-1], collapse = "+")
  log_model <- glm(paste0("as.factor(Group)~", formula),
                   data = data1, family = binomial(link = "logit"))
  glm.probs <- predict(log_model, data1[, -1, drop = F], type = "response")
  
  data1[, "Expression"] <- glm.probs
  data1[, "Feature"] <- metaname
  data1 <- data1[, c("Feature", "Group", "Expression")]
  data2 <- data1[, c("Expression"), drop = F]
  colnames(data2) <- metaname
  data2 <- as.data.frame(t(data2))
  
  auc <- auto_plotroc(data = data1,
                      mapname = mapname,
                      savepath = savepath,
                      ...)
  
  auc$auc[, "ROC图名称"] <- mapname
  
  auc[["data"]] <- data2
  
  # print("auto_logroc运行完成")
  if (saveroc) {
    savexlsx3(data = auc$auc, 
              filename = paste0(savepath,"/AUC.xlsx"), 
              sheet = mapname)
    savexlsx3(data = data2, 
              filename = paste0(savepath,"/Logistic.xlsx"), 
              sheet = mapname)
  } else {
    return(auc)
  }
}

#' 逻辑回归后roc可视化
#'
#' @param data 数据
#' @param name 保存文件名
#' @param type 保存图片格式
#' @param saveroc 保存roc数据
#' @param num 图上绘制roc数据
#' @param all 逻辑，是否一起绘制
#' @param ... 见[auto_logroc()]
#'
#' @export
auto_logrocnum <- function(data,
                           mapname = "LogsticROC",
                           num = dim(data)[1]-1,
                           all = F,
                           saveroc = T,
                           savepath ="./",
                           ...) {
  
  if (row.names(data)[1] != "Group") {
    stop("请在数据第二行添加Group分组信息")
  }
  
  if(num == 0){
    num <- dim(data)[1]-1
  }
  
  rownumber <- 2:nrow(data)
  auc <- NULL
  data2 <- NULL
  
  if (length(rownumber) == 1 | num == 1) {
    stop("选择物质数量少于2，请最少选择2个的物质")
  }
  
  for (i in 1:length(num)) {
    if (length(rownumber) < num[i]) {
      # print(num[i])
      stop("组合特征数量少于需组合数量")
    }
    combndata <- combn(rownumber, num[i])
    
    for (j in 1:ncol(combndata)) {
      aucdata <- auto_logroc(data = data[c(1, combndata[, j]), ],
                             mapname = paste0(mapname, "-", num[i], "-", j),
                             saveroc = F,
                             savepath = savepath,
                             ...)
      
      if (all & length(combndata[, j]) > 1) {
        data3 <- rbind(data[1, ], aucdata$data, data[combndata[, j], ])
        
        aucdata2 <- auto_roc(data = data3,
                             mapname = paste0(mapname, "-", num[i], "-", j, "-all"),
                             number = NULL,
                             saveroc = F,
                             savepath = savepath,
                             ...)
        
        aucdata$plot <- aucdata2$plot
        auc <- rbind(auc, aucdata2$auc)
      }
      
      auc <- rbind(auc, aucdata$auc)
      data2 <- rbind(data2, aucdata$data)
    }
  }
  
  aucdata$auc <- auc
  aucdata$data <- data2
  
  # print("auto_logrocnum运行完成")
  
  if (saveroc) {
    savexlsx1(data = auc, 
              filename = paste0(savepath,"/AUC.xlsx"), 
              sheet = mapname)
    savexlsx3(data = data2, 
              filename = paste0(savepath,"/Logistic.xlsx"), 
              sheet = mapname)
  } else {
    return(aucdata)
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_common_logroc <- map_autodraw$new(moudle = auto_logrocnum,row.names = 1)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "表达数据矩阵",
                      required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")

  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "分析结果", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("jpg","pdf","html"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf","html"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 特有参数
  parser$add_argument("-nsr","--nosaveroc", default = T, action = "store_false", 
                      help="是否保存roc结果数据",dest = "saveroc")
  parser$add_argument("-na","--nAUC", default = T,  action = "store_false", 
                      help="是否显示AUC",dest = "AUC")
  parser$add_argument("-n","--num", default = 0, type = "integer", 
                      help="随机组合进行逻辑回归roc计算")
  parser$add_argument("-al","--all", default = F, action = "store_true", 
                      help="是否一起绘制随机组合进行逻辑回归的roc")
  parser$add_argument("-z","--zip",default = F, help = "是否压缩",action='store_true')
  parser$add_argument("-s","--savepath",default = "分析结果", help = "结果输出路径")
  
  args <- parser$parse_args()
  
  zip <- args$zip
  args$zip <- NULL
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}

  rocresult <- do.call(what = map_common_logroc, args = args)
  
  if(zip){
    zip::zip(zipfile = "分析结果.zip",files = args$saveptah)
    unlink(list.files(pattern = "[^z][^i][^p]$",include.dirs = T,all.files = T), recursive = T)
  }
  
}

#' 根据文件进行logroc可视化
#' 
#' @export
map_common_logroc <- map_autodraw$new(moudle = auto_logrocnum,row.names = 1)$draw
