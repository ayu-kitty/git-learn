#!/opt/conda/bin/Rscript

#' 时间序列/趋势分析
#'
#' @param data 数据
#' @param name 保存文件名
#' @param circle_num 循环次数
#' @param max_cluster 最大聚类数
#' @param cluster 指定聚类数
#' @param filter_std 逻辑，是否筛选
#' @param standardise 逻辑，是否标准化
#' @param ... 见[auto_stemanalyst()]
#'
#' @export
auto_stem <- function(filename,
                      resultPath = "时间序列及趋势分析",
                      circle_num = 5,
                      max_cluster = 50,
                      cluster = NULL,
                      filter_std = FALSE,
                      standardise = T,
                      ...) {
  
  wd <- getwd()
  
  data <- readdata(filename = filename, row.names = 1)
  
  setwddir(filename = resultPath,force = T)
  
  auto_stemanalyst(data = data,
                   circle_num = circle_num,
                   max_cluster = max_cluster,
                   cluster = cluster,
                   filter_std = filter_std,
                   standardise = standardise,
                   ...)
  
  # 生成报告
  mkreport(type = "时间序列及趋势分析")
  
  setwd(wd)
  print("STEM分析结束，热图、聚类图保存结束")
}

# 设置循环次数，获得最佳聚类个数
#' @export
cir_num <- function(circle_num = 2,
                    max_cluster = 50,
                    eset = eset,
                    m = m) {
  print(paste("circle_num为：", circle_num, " ; max_cluster为：", max_cluster, sep = ""))
  if (circle_num > 6) {
    circle_num <- 6
  } else {
    circle_num <- circle_num
  }
  c_list <- list()
  
  pdf("Dmin.pdf")
  for (j in seq(1, circle_num)) {
    tmp <- Mfuzz::Dmin(eset,
                       m = m,
                       crange = seq(10, max_cluster, 5),
                       repeats = 3,
                       visu = TRUE)
    
    for (i in seq(1, (length(tmp) - 1))) {
      if ((tmp[i] - tmp[i + 1]) < 0.005) {
        c1 <- 5 + 5 * i
        break
      } else {
        c1 <- max_cluster
      }
    }
    c_list <- append(c_list, c1)
  }
  c <- round(rowMeans(as.data.frame(c_list)))
  dev.off()
  return(c)
}

# 画cluster的折线图
#' @export
cl_draw <- function(data,
                    choosedata,
                    i,
                    standardise,
                    width = 9,
                    height = 7,
                    ...) {
  suppressMessages(library("ggplot2"))
  theme_set(theme_classic())
  p <- ggplot(data = choosedata, aes(x = factor(choosedata$Samples, levels = colnames(data)),
                                     y = value)) +
    geom_line(aes(color = Membership, group = Membership)) +
    labs(title = paste0("Cluster ", i),
         subtitle = paste("contains ", length(choosedata$Membership) / length(colnames(data)), " Feature" , sep = ""),
         x = "Samples") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.text.align = 0.5)
  
  col <- c("#FFFF66", "#FFFF33", "#FFFF00", "#33FF66", "#33CC33", "#33CC00",
           "#0099FF", "#0066FF", "#6633CC", "#663399", "#663366", "#990099",
           "#FF0066", "#FF0033", "#CC0000")
  
  if (!standardise) {
    p <- p + scale_color_gradientn(colours = col, breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(name = "Expression Change (log2(v(i)/v(0))")
  } else {
    p <- p + scale_color_gradientn(colours = col, breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(name = "Expression Abundance (Z-score)")
  }
  
  # 保存图片
  print(paste("保存Cluster ", i, "的聚类图", sep = ""))
  ggplotsave(plot = p,
             mapname = paste0("Cluster ", i),
             width = width, height = height,
             ...)
}

#' @export
auto_stemanalyst <- function(data,
                             circle_num = 5,
                             max_cluster = 50,
                             cluster = NULL,
                             filter_std = FALSE,
                             min.std = 0,
                             standardise = T,
                             ...) {
  print("数据预处理")
  data1 <- data
  colnames(data1) <- paste0(colnames(data1), ".raw")
  colone_value <- data1[, 1]
  if (!standardise) {
    print("数据进行log2g归一化处理")
    
    value_log <- function(x, na.rm = FALSE) (log2(x / colone_value))
    data2 <- dplyr::mutate_all(.tbl = data1, .funs = value_log)
    colnames(data2) <- paste0(colnames(data), ".Expression Change(log2(v(i)/v(0))")
    row.names(data2) <- rownames(data1)
  } else {
    print("数据进行Z-score标准化处理...")
    data2 <- t(scale(t(data1))) # 对row进行Z-score标准化
    colnames(data2) <- paste0(colnames(data), ".Expression Abundance(Z-score)")
  }
  
  print("热图总图保存")
  data3 <- as.data.frame(data2)
  names(data3) <- colnames(data)
  auto_heatmap(data3, 
               mapname = "Summary-heatmap",
               ...)
  
  print("数据保存")
  
  data4 <- cbind(Feature = row.names(data2),
                 data2,
                 data1)
  
  savexlsx1(data = data4,
            filename = "分析原始数据.xlsx",
            sheet = "分析原始数据")
  
  # 需要matrix数据类型
  count_matrix <- data.matrix(data2)
  
  suppressMessages(library("Mfuzz"))
  suppressMessages(library("dplyr"))
  eset <- new("ExpressionSet", exprs = count_matrix)
  
  # 根据标准差去除样本间差异太小的基因
  if (!filter_std) {
    eset <- eset
  } else {
    print("当前数据正在进行过滤处理...")
    eset <- Mfuzz::filter.std(eset, min.std = min.std)
  }
  
  print("开始分析数据")
  #  评估出最佳的m值
  m <- Mfuzz::mestimate(eset)
  
  # 判断最佳的聚类个数
  if (max_cluster > nrow(data1)) {
    warning("max_cluster数大于矩阵行数", immediate. = T)
    max_cluster <- nrow(data1)
  }
  
  if (!is.null(cluster)) {
    c <- cluster
    print(paste0("聚类个数修改为: ", c))
  } else {
    c <- cir_num(circle_num = circle_num,
                 m = m,
                 eset = eset,
                 max_cluster = max_cluster,
                 ...)
    print(paste0("获得最佳聚类个数为: ", c))
  }
  
  # 聚类总图
  n1 <- ceiling(sqrt(c + 1))
  n2 <- ceiling((c + 1) / n1)
  
  cl <- Mfuzz::mfuzz(eset, c = c, m = m)
  
  plotfile(mapname = "Cluster-summary",
           height = n1 * 4,
           width = n2 * 6,
           ...)
  Mfuzz::mfuzz.plot(eset, cl,
                    mfrow = c(n2, n1),
                    new.window = FALSE,
                    time.labels = names(data))
  
  mfuzzColorBar <- function (col, horizontal = FALSE, ...){
    require(marray) || stop("Library marray is required")
    if (missing(col)) {
      col <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700", 
                        "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00", 
                        "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", 
                        "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", 
                        "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7", 
                        "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", 
                        "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", 
                        "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF", 
                        "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF", 
                        "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", 
                        "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", 
                        "#FF0060", "#FF0048", "#FF0030", "#FF0018")
    }else if (length(col) > 1) {
    }else if (col == "fancy") {
      fancy.blue <- c(c(255:0), rep(0, length(c(255:0))), 
                      rep(0, length(c(255:150))))
      fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
      fancy.red <- c(c(0:255), rep(255, length(c(255:0))), 
                     c(255:150))
      col <- rgb(b = fancy.blue/255, g = fancy.green/255, 
                 r = fancy.red/255)
    }
    par(mar = c(5, 2, 4, 3) + 0.1)
    maColorBar(seq(0, 1, 0.01), col = col, horizontal = FALSE, 
               k = 11, ...)
  }
  
  suppressWarnings(mfuzzColorBar(main = "Membership", cex.main = 1)) # 画lengend，默认dpi为300
  plotsave()
  
  # 保存每个cluster中的element个数
  
  cl.size <- data.frame("Cluster number" = paste0("Cluster ", seq(1, c)),
                        "Cluster amount" = cl$size, check.names = F)
  savexlsx1(data = cl.size,
            filename = "Cluster Statistics.xlsx",
            sheet = "Cluster Statistics")
  
  
  # 查看基因和cluster之间的membership
  membership <- cl$membership
  membership1 <- data.frame(membership)
  names(membership1) <- paste0("Cluster ", seq(1, c))
  
  membership1 <- cbind(Feature = row.names(membership1),
                       membership1)
  
  savexlsx1(data = membership1,
            filename = "Membership.xlsx",
            sheet = "Membership")
  
  # 提取cluster下的element,并保存成文件
  for (i in seq(1, c)) {
    cluster_i <- data.frame("Cluster number" = cl$cluster[cl$cluster == i],
                            check.names = F)
    cluster_i[, 1] <- paste0("Cluster ", cluster_i[, 1])
    
    print(paste0("保存Cluster ", i, "的热图"))
    
    data_i <- data3[row.names(cluster_i), , drop = F] # 筛选出该Cluster的metabolites
    # 画热图
    auto_heatmap(data = data_i,
                 mapname = paste0("Cluster ", i, "-heatmap"),
                 ...)
    
    cluster_i <- cbind(cluster_i, data4[row.names(cluster_i), ]) # 合并
    # 提取每个cluster的element的membership
    mem.ship <- data.frame("Membership" = membership[row.names(cluster_i), i])
    cluster_i <- cbind(cluster_i, mem.ship) # 合并
    # 保存分析结果
    print(paste0("保存Cluster ", i, "的结果"))
    savexlsx1(data = cluster_i,
              filename = paste0("Cluster ", i, ".xlsx"),
              sheet = paste0("Cluster ", i))
    
    # 查看属于cluster cores的基因list. cluster cores为membership > 0.7的metabolites
    acore.list <- Mfuzz::acore(eset, cl, min.acore = 0.7)
    acore.list[[i]] <- dplyr::rename(acore.list[[i]], Feature = NAME)
    acore.list[[i]] <- dplyr::rename(acore.list[[i]], Membership = MEM.SHIP)
    
    savexlsx1(data = acore.list[[i]],
              filename = paste0("Cluster ", i, ".xlsx"),
              sheet = paste0("Cluster ", i, " cores"))
    
    # 画图，线的颜色以各个代谢物在该cluster的membership来定
    
    choosedata <- data.frame(t(select(cluster_i, contains("Expression")))) # 抽提画图信息（时间轴，表达量）
    h <- length(rownames(choosedata))
    choosedata <- reshape2::melt(choosedata) # 将宽数据框变为长数据框， 便于画图
    # choosedata['Expression Change(log2(v(i)/v(0))'] <- choosedata['value']
    choosedata["Samples"] <- colnames(data)
    choosedata["Membership"] <- rep(cluster_i$Membership, each = h, times = 1) # 添加membership，指示颜色
    
    cl_draw(data = data,
            choosedata = choosedata,
            i = i,
            standardise = standardise,
            ...)
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-f","--filename",default = "表达数据矩阵", nargs = "+",
                      help = "表达数据矩阵文件")
  parser$add_argument("-cn","--circle_num", default = 5, type = "integer", help="循环次数")
  parser$add_argument("-mc","--max_cluster", default = 30, type = "integer", help="最大聚类数")
  parser$add_argument("-c","--cluster", default = 10, type = "integer", help="指定聚类数")
  parser$add_argument("-fs","--filter_std", default = F, action = "store_true", help="逻辑，是否筛选")
  parser$add_argument("-nst","--nstandardise", default = T,  action = "store_false", help="逻辑，是否标准化",
                      dest = "standardise")
  
  parser$add_argument("-r","--resultPath",default = "时间序列分析", help = "结果输出路径")
  parser$add_argument("-z","--zip",default = F, help = "是否压缩",action='store_true')
  
  args <- parser$parse_args()
  
  zip <- args$zip
  args$zip <- NULL
  
  stemresult <- do.call(what = auto_stem, args = args)
  
  if(zip){
    zip::zip(zipfile = "分析结果.zip",files = args$resultPath)
    unlink(list.files(pattern = "[^z][^i][^p]$",include.dirs = T,all.files = T), recursive = T)
  }
}
