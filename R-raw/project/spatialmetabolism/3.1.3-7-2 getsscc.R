#!/opt/conda/bin/Rscript

#' 获取sscc分析特征值及将特征值进行可视化
#'
#' @param datapath sscc数据路径
#' @param samplename 样本名称
#' @param savename 保存名称
#' @param mode 正负离子模式
#' @param qualitativepath 定性数据路径
#' @param savepath 保存路径
#' @param moudle 可视化模块
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param ... 见moudle函数
#'
#' @export
getssccfeature <- function(datapath = "./sample/cluster/sscc/",
                           samplename = NULL,
                           savename = samplename,
                           mode = "neg",
                           qualitativepath = "./sample/qualitative/Qualitative.xlsx",
                           savepath = "./",
                           moudle = NULL,
                           mapname = "",
                           imagetype = c("png", "pdf"),
                           filename = paste0(datapath, samplename, "-", mode, "-data.rds"),
                           ...) {
  suppressMessages(library("Cardinal"))
  
  # rds <- paste0(datapath, samplename, "-", mode, "-data.rds")
  
  data <- readdata(filename = filename)
  
  ssc <- data$sscc
  
  if(file.exists(qualitativepath)){
    qdata <- readdata(filename = qualitativepath,
                      sheet = mode)
    
    if (is.numeric(qdata$mz)) {
      qdata$mz <- format(qdata$mz, nsmall = 5, trim = T)
    }
  }else{
    qdata <- NULL
  }

  cordata <- ssc@resultData@listData[[1]][["statistic"]]
  colnames(cordata) <- paste0("Cluster ", colnames(cordata))
  rownames(cordata) <- format(ssc@featureData@mz, nsmall = 5, trim = T)
  
  ssc_top <- as.data.frame(cordata)
  ssc_top[, "mz"] <- format(ssc@featureData@mz, nsmall = 5, trim = T)
  if(!is.null(qdata)){
    ssc_top <- merge(x = ssc_top, y = qdata, by.x = "mz", by.y = "mz", all.x = T)
  }
  
  savexlsx1(data = ssc_top,
            filename = paste0(savepath, "Cluster-特征值.xlsx"),
            sheet = paste0(savename, "-", mode))
  
  if (is.function(moudle)) {
    returndata <- moudle(cordata,
                         mapname = paste0(mapname, savename, "-", mode),
                         savepath = savepath,
                         imagetype = imagetype,
                         ...)
  }
}


#' 获取sscc分析特征值及将特征值进行相关性可视化
#'
#' @param moudle 可视化模块
#' @param mapname 图片名称
#' @param trans 数据是否转置
#' @param order 相关性图是否排序
#' @param title 标题
#' @param ... 见moudle函数
#'
#' @export
getssccfeaturetocor <- function(moudle = datatocorrplot,
                                trans = F,
                                order = "original",
                                title = "",
                                mapname = "Cluster-correlation-",
                                annotation_colors = NULL,
                                ...) {
  
  getssccfeature(moudle = moudle,
                 transx = trans,
                 order = order,
                 title = title,
                 mapname = mapname,
                 ...)
  
}


#' 获取sscc分析的特征物质及将特征物质进行可视化
#'
#' @param datapath sscc数据路径
#' @param samplename 样本名称
#' @param savename 保存名称
#' @param mode 正负离子模式
#' @param qualitativepath 定性数据路径
#' @param savepath 保存路径
#' @param moudle 可视化模块
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param topFeaturesn 数量，前几的特征物质
#' @param savetopdata 逻辑，是否保存特征物质数据
#' @param annosamplename 逻辑,是否注释样本名
#' @param sampletogroup 逻辑，是否sample改成Group
#' @param annocluster 逻辑，是否注释聚类组
#' @param classtogroup 逻辑，是否class改成Group
#' @param annometabo 逻辑，是否注释代谢物名
#' @param order 逻辑，是否根据数量排序
#' @param annotation_colors 颜色
#' @param ... 见moudle函数
#'
#' @export
getsscctopfeature <- function(datapath = "./sample/cluster/sscc/",
                              samplename = NULL,
                              savename = samplename,
                              mode = "neg",
                              topFeaturesn = 10,
                              topFeaturesclustersave = F,
                              savetopdata = T,
                              qualitativepath = "./sample/qualitative/Qualitative.xlsx",
                              savepath = "./",
                              moudle = NULL,
                              mapname = "",
                              imagetype = c("png", "pdf"),
                              annosamplename = T,
                              sampletogroup = F,
                              annocluster = T,
                              classtogroup = F,
                              annometabo = F,
                              order = "cor",
                              annotation_colors = list(Cluster = SelectColors(palette = "customecol2", n = 50),
                                                       Sample = SelectColors(palette = "cold", n = 30)),
                              filename = paste0(datapath, samplename, "-", mode, "-data.rds"),
                              ...) {
  suppressMessages(library("Cardinal"))
  # rds <- paste0(datapath, samplename, "-", mode, "-data.rds")
  data <- readdata(filename = filename)
  
  spectradata <- data$spectradata
  ssc <- data$sscc
  
  if(file.exists(qualitativepath)){
    
    qdata <- readdata(filename = qualitativepath,sheet = mode)
    
    if (is.numeric(qdata$mz)) {
      qdata$mz <- format(qdata$mz, nsmall = 5, trim = T)
    }
    
    if (paste0(mode, "-all") %in% getsheetname(qualitativepath)) {
      qdata_all <- readdata(filename = qualitativepath,sheet = paste0(mode, "-all"))
      
      if (is.numeric(qdata_all$mz)) {
        qdata_all$mz <- format(qdata_all$mz, nsmall = 5, trim = T)
      }
      
    } else {
      qdata_all <- NULL
    }
    
  }else{
    qdata <- NULL
    qdata_all <- NULL
  }
  
 
  
  createdir(filename = savepath)
  
  if (!is.null(topFeaturesn)) {
    num <- max(as.numeric(ssc@resultData@listData[[1]][["class"]]))
    ssc_top <- NULL
    for (i in 1:num) {
      ssc_top2 <- as.data.frame(topFeatures(ssc, class == i, n = topFeaturesn))
      ssc_top2 <- ssc_top2[ssc_top2$statistic > 0,,drop = F]
      ssc_top <- rbind(ssc_top, ssc_top2)
    }
    
    ssc_top$mz <- format(ssc_top$mz, nsmall = 5, trim = T)
    
    ssc_top_all <- ssc_top
    if(!is.null(qdata)){
      if (is.numeric(qdata$mz)) {
        qdata$mz <- format(qdata$mz, nsmall = 5, trim = T)
      }
      ssc_top <- merge(x = ssc_top, y = qdata, by.x = "mz", by.y = "mz", all.x = T)
    }
    ssc_top <- ssc_top[order(abs(ssc_top$statistic), decreasing = T), ]
    ssc_top <- ssc_top[order(as.numeric(ssc_top$class)), ]
    
    if (savetopdata) {
      
      savexlsx1(data = ssc_top,
                filename = paste0(savepath, "Cluster-特征物质-top", topFeaturesn, ".xlsx"),
                sheet = paste0(savename, "-", mode))
      
      if (topFeaturesclustersave) {
        for (clustenum in unique(ssc_top$class)) {
          
          savexlsx1(data = ssc_top[ssc_top$class == clustenum, ],
                    filename = paste0(savepath, "Cluster-特征物质-top", topFeaturesn, "-cluster.xlsx"),
                    sheet = paste0(savename, "-", mode, "-", clustenum))
          
        }
      }
      
      if (!is.null(qdata_all)) {
        
        ssc_top_all <- merge(x = ssc_top_all, y = qdata_all, by.x = "mz", by.y = "mz", all.x = T)
        ssc_top_all <- ssc_top_all[order(abs(ssc_top_all$statistic), decreasing = T), ]
        ssc_top_all <- ssc_top_all[order(as.numeric(ssc_top_all$class)), ]
        
        savexlsx1(data = ssc_top_all,
                  filename = paste0(savepath, "Cluster-特征物质-top", topFeaturesn, ".xlsx"),
                  sheet = paste0(savename, "-", mode, "-all"))
        
        if (topFeaturesclustersave) {
          for (clustenum in unique(ssc_top_all$class)) {
            
            savexlsx1(data = ssc_top_all[ssc_top_all$class == clustenum, ],
                      filename = paste0(savepath, "Cluster-特征物质-top", topFeaturesn, "-cluster-all.xlsx"),
                      sheet = paste0(savename, "-", mode, "-", clustenum))
            
          }
        }
      }
    }
    
    ssc_top <- ssc_top[order(abs(ssc_top$statistic), decreasing = T), ]
    ssc_top <- ssc_top[!duplicated(ssc_top$mz), ]
    
    if (order == "num") {
      orderlevel <- table(ssc@resultData@listData[[1]][["class"]])
      orderlevel <- orderlevel[order(orderlevel, decreasing = T)]
      levelclass <- factor(ssc_top$class, level = names(orderlevel))
      ssc_top <- ssc_top[order(levelclass), ]
    }else if (order == "cor") {
      orderlevel <- table(ssc@resultData@listData[[1]][["class"]])
      orderlevel <- orderlevel[order(orderlevel, decreasing = T)]
      
      cordata <- ssc@resultData@listData[[1]][["statistic"]]
      cordata <-  as.data.frame(cordata)
      cordata2 <- corrcal(x = as.data.frame(cordata),transx = F)
      linkdata <- cordata2[["framedata"]][["linkdata"]]
      
      for ( k in 2:(length(orderlevel)-1)) {
        linkdata2 <- linkdata[linkdata$Featurex == names(orderlevel)[k-1] | linkdata$Featurey == names(orderlevel)[k-1],]
        linkdata <- linkdata[!(linkdata$Featurex == names(orderlevel)[k-1] | linkdata$Featurey == names(orderlevel)[k-1]),]
        linkdata2 <- linkdata2[which.max(linkdata2$cor),]
        linkdata2 <- t(linkdata2)[1:2,1]
        linkdata2 <- linkdata2[ linkdata2!= names(orderlevel)[k-1]]
        orderlevel <- c(orderlevel[1:(k-1)],orderlevel[which(names(orderlevel) == linkdata2)],orderlevel[-c(1:(k-1),which(names(orderlevel) == linkdata2))])
      }
 
      levelclass <- factor(ssc_top$class, level = names(orderlevel))
      ssc_top <- ssc_top[order(levelclass), ]
    }else {
      ssc_top <- ssc_top[order(as.numeric(ssc_top$class)), ]
    }
    
    spectradata2 <- spectradata[as.character(ssc_top$mz), ]
  } else {
    spectradata2 <- spectradata
  }
  
  if (is.function(moudle)) {
    
    if (annometabo) {
      annometaboname <- whto(qdata[, c("mz", "Metabolites")], row.names(spectradata2))
      row.names(spectradata2)[!is.na(annometaboname)] <- annometaboname[!is.na(annometaboname)]
    }
    spectradata3 <- spectradata2[0, ]
    
    annotation_colors2 <- list()
    
    if (annosamplename & length(unique(data$plotdata$Sample)) > 1) {
      samplecolor <- annotation_colors$Sample[1:length(unique(data$plotdata$Sample))]
      names(samplecolor) <- unique(data$plotdata$Sample)
      
      if (sampletogroup) {
        spectradata3["Group", ] <- data$plotdata$Sample
        annotation_colors2[["Group"]] <- samplecolor
      } else {
        spectradata3["Sample", ] <- data$plotdata$Sample
        annotation_colors2[["Sample"]] <- samplecolor
      }
    }
    
    if (annocluster) {
      
      clustercolor <- annotation_colors$Cluster[1:length(levels(ssc@resultData@listData[[1]][["class"]]))]
      names(clustercolor) <- levels(ssc@resultData@listData[[1]][["class"]])
      spectradata3 <- as.data.frame(spectradata3)
      
      if (classtogroup) {
        spectradata3["Group", ] <- ssc@resultData@listData[[1]][["class"]]
        annotation_colors2[["Group"]] <- clustercolor
      } else {
        spectradata3["Cluster", ] <- ssc@resultData@listData[[1]][["class"]]
        annotation_colors2[["Cluster"]] <- clustercolor
      }
      print("2")
      spectradata2 <- rbind(spectradata3, spectradata2)
      # row.names(spectradata2) <- c("Cluster", format(ssc_top$mz, nsmall = 5,trim = T))
      if (order == "num") {
        orderlevel <- table(ssc@resultData@listData[[1]][["class"]])
        orderlevel <- orderlevel[order(orderlevel, decreasing = T)]
        levelclass <- factor(ssc@resultData@listData[[1]][["class"]], level = names(orderlevel))
        spectradata2 <- spectradata2[, order(levelclass)]
      }else if (order == "cor") {
        orderlevel <- table(ssc@resultData@listData[[1]][["class"]])
        orderlevel <- orderlevel[order(orderlevel, decreasing = T)]
        
        cordata <- ssc@resultData@listData[[1]][["statistic"]]
        cordata <-  as.data.frame(cordata)
        cordata2 <- corrcal(x = as.data.frame(cordata),transx = F)
        linkdata <- cordata2[["framedata"]][["linkdata"]]
        
        for ( k in 2:(length(orderlevel)-1)) {
          linkdata2 <- linkdata[linkdata$Featurex == names(orderlevel)[k-1] | linkdata$Featurey == names(orderlevel)[k-1],]
          linkdata <- linkdata[!(linkdata$Featurex == names(orderlevel)[k-1] | linkdata$Featurey == names(orderlevel)[k-1]),]
          linkdata2 <- linkdata2[which.max(linkdata2$cor),]
          linkdata2 <- t(linkdata2)[1:2,1]
          linkdata2 <- linkdata2[ linkdata2!= names(orderlevel)[k-1]]
          orderlevel <- c(orderlevel[1:(k-1)],orderlevel[which(names(orderlevel) == linkdata2)],orderlevel[-c(1:(k-1),which(names(orderlevel) == linkdata2))])
        }
        
        levelclass <- factor(ssc@resultData@listData[[1]][["class"]], level = names(orderlevel))
        spectradata2 <- spectradata2[, order(levelclass)]
      }else {
        spectradata2 <- spectradata2[, order(ssc@resultData@listData[[1]][["class"]])]
      }
    } else {
      spectradata2 <- rbind(spectradata3, spectradata2)
    }
    
    returndata <- moudle(spectradata2,
                         mapname = paste0(mapname, savename, "-", mode),
                         savepath = savepath,
                         imagetype = imagetype,
                         annotation_colors = annotation_colors2,
                         ...)
    
    write.csv(x = spectradata2,
              file = paste0(savepath, mapname, savename, "-", mode, ".csv"))
    
  }
}

#' 获取sscc分析的特征物质及将特征物质进行热图可视化
#'
#' @param moudle 可视化模块
#' @param mapname 图片名称
#' @param annosamplename 逻辑,是否注释样本名
#' @param annocluster 逻辑，是否注释聚类组
#' @param annotation_colors 注释颜色
#' @param colgroup 是否注释
#' @param cluster_rows 是否横聚类
#' @param cluster_cols 是否纵聚类
#' @param show_colnames 是否显示样本名称
#' @param breaks 图例范围
#' @param color 热图颜色
#' @param ... 见moudle函数
#'
#' @export
getsscctopfeaturetoheatmap <- function(moudle = auto_heatmap,
                                       annosamplename = F,
                                       annocluster = T,
                                       colgroup = ifelse(annocluster & annosamplename, 2, ifelse(annocluster | annosamplename, 1, 0)),
                                       cluster_rows = F,
                                       cluster_cols = F,
                                       show_colnames = F,
                                       breaks = c(-500:500) / 100,
                                       color = colorRampPalette(c("violetred", "black", "gold"))(1001),
                                       mapname = "Cluster-heatmap-",
                                       annotation_colors = list(Cluster = SelectColors(palette = "customecol2", n = 50),
                                                                Sample = SelectColors(palette = "cold", n = 30)),
                                       ...) {
  
  getsscctopfeature(moudle = moudle,
                    annosamplename = annosamplename,
                    annocluster = annocluster,
                    colgroup = colgroup,
                    cluster_rows = cluster_rows,
                    cluster_cols = cluster_cols,
                    show_colnames = show_colnames,
                    breaks = breaks,
                    color = color,
                    mapname = mapname,
                    annotation_colors = annotation_colors,
                    ...)
  
}


#' 获取sscc分析的特征物质及将特征物质进行小提琴可视化
#'
#' @param moudle 可视化模块
#' @param mapname 图片名称
#' @param savetopdata 逻辑，是否保存特征物质数据
#' @param annosamplename 逻辑,是否注释样本名
#' @param annocluster 逻辑，是否注释聚类组
#' @param classtogroup 逻辑，是否class改成Group
#' @param order 逻辑，是否根据数量排序
#' @param annotation_colors 颜色
#' @param boxmoudle 小提琴模块
#' @param aspect.ratio 图像长宽比
#' @param height 图片高度
#' @param width 图片宽度
#' @param legend 是否显示图例
#' @param xlab x轴名称
#' @param ... 见moudle函数
#'
#' @export
getsscctopfeaturetoviolin <- function(moudle = auto_boxchart,
                                      boxmoudle = auto_violin2,
                                      annosamplename = F,
                                      annocluster = annocluster,
                                      classtogroup = T,
                                      aspect.ratio = 9 / 16,
                                      height = 5,
                                      width = 9,
                                      legend = "none",
                                      xlab = "Cluster",
                                      mapname = "Cluster-Violin-all-",
                                      order = "none",
                                      savetopdata = F,
                                      annotation_colors = list(Cluster = SelectColors(palette = "customecol2", n = 50),
                                                               Sample = SelectColors(palette = "cold", n = 30)),
                                      other = list(theme_classic(),
                                                   theme(strip.background = element_rect(color = "white", fill = "white"),
                                                         panel.grid = element_blank(),
                                                         plot.margin = unit(c(1,1,1,1),"cm"),
                                                         legend.position = "none")),
                                      ...) {
  
  getsscctopfeature(moudle = moudle,
                    boxmoudle = boxmoudle,
                    annosamplename = annosamplename,
                    classtogroup = classtogroup,
                    mapname = mapname,
                    aspect.ratio = aspect.ratio,
                    height = height,
                    width = width,
                    legend = legend,
                    xlab = xlab,
                    order = order,
                    savetopdata = savetopdata,
                    annotation_colors = annotation_colors,
                    other = other,
                    ...)
  
}

#' 获取sscc的成像图
#'
#' @param datapath sscc数据路径
#' @param samplename 样本名称
#' @param savename 保存名称
#' @param mode 正负离子模式
#' @param savepath 保存路径
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param annotation_colors 颜色
#' @param lightmode 成像模式
#' @param lengedallshow 逻辑，聚类图例全部显示
#' @param ... 见[clusterlimage()]
#'
#' @export
getsscctoimage <- function(datapath = "./sample/cluster/sscc/",
                           samplename = NULL,
                           savename = samplename,
                           mode = "neg",
                           savepath = "./",
                           mapname = "Cluster-sscc-",
                           imagetype = c("png", "pdf"),
                           lightmode = T,
                           lengedallshow = T,
                           annotation_colors = list(Cluster = SelectColors(palette = "customecol2", n = 50),
                                                    Sample = SelectColors(palette = "cold", n = 30)),
                           filename = paste0(datapath, samplename, "-", mode, "-data.rds"),
                           ...) {
  
  if (!file.exists(filename)) {
    print(paste0(filename, "文件不存在"))
    return()
  }
  
  imzmlclusterimage(filename = filename,
                    savepath = savepath,
                    mapname = paste0(mapname, savename, "-", mode),
                    imagetype = imagetype,
                    lightmode = lightmode,
                    lengedallshow = lengedallshow,
                    col = annotation_colors$Cluster,
                    ...)
  
}

#' 获取sscc的质谱图
#'
#' @param datapath sscc数据路径
#' @param samplename 样本名称
#' @param savename 保存名称
#' @param mode 正负离子模式
#' @param savepath 保存路径
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param annotation_colors 颜色
#' @param lightmode 成像模式
#' @param ... 见[clusterlplot()]
#'
#' @export
getsscctoplot <- function(datapath = "./sample/cluster/sscc/",
                          samplename = NULL,
                          savename = samplename,
                          mode = "neg",
                          savepath = "./",
                          mapname = "Cluster-tstatistics-",
                          imagetype = c("png", "pdf"),
                          lightmode = T,
                          annotation_colors = list(Cluster = SelectColors(palette = "customecol2", n = 50),
                                                   Sample = SelectColors(palette = "cold", n = 30)),
                          filename = paste0(datapath, samplename, "-", mode, "-data.rds"),
                          ...) {

  if (!file.exists(filename)) {
    print(paste0(filename, "文件不存在"))
    return()
  }
  
  imzmlclusterplot(filename = filename,
                   savepath = savepath,
                   mapname = paste0(mapname, savename, "-", mode),
                   imagetype = imagetype,
                   lightmode = lightmode,
                   col = annotation_colors$Cluster,
                   ...)
}
