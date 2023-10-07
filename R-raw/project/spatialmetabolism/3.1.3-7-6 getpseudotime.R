#!/opt/conda/bin/Rscript

#' getpseudotime_data
#'
#' 获取拟时序数据
#'
#' @param datapath 拟时序数据路径
#' @param samplename 样本名称
#' @param savename 保存名称
#' @param mode 正负离子模式
#' @param savepath 保存路径
#' @param name 保存文件名
#'
#' @export
getpseudotime_data <- function(datapath = "./sample/cluster/pseudotime/",
                               savepath = "./",
                               mode = "neg",
                               filename = paste0(datapath,"/",samplename, "-", mode, "-data.rds"),
                               samplename = "1462",
                               savename = samplename,
                               name = "PseudotimeData.xlsx") {
  # 加载monocle
  suppressMessages(library("monocle"))
  
  # rds <- paste0(datapath,"/",samplename, "-", mode, "-data.rds")
  
  data <- readdata(filename = filename)
  
  pdata <- pData(data$cds)

  savexlsx1(data = pdata,
            filename = paste0(savepath,"/",name),
            sheet = mode)
}

#' getpseudotime_plotcell
#'
#' 获取拟时序图
#'
#' @param datapath 拟时序数据路径
#' @param samplename 样本名称
#' @param savename 保存名称
#' @param mode 正负离子模式
#' @param mapname 图片名称
#' @param filtermz 绘制mz
#' @param mapmoudle 绘图函数
#' @param ... 见mapmoudle
#'
#' @export
getpseudotime_plotcell <- function(datapath = "./sample/cluster/pseudotime/",
                                   samplename = "1462",
                                   savename = samplename,
                                   mode = "neg",
                                   filename = paste0(datapath,"/",samplename, "-", mode, "-data.rds"),
                                   filtermz = NULL,
                                   mapname = "",
                                   mapmoudle = getpseudotime_plotcelltrajectory ,
                                   ...) {
  # 加载monocle
  suppressMessages(library("monocle"))

  # rds <- paste0(datapath,"/",samplename, "-", mode, "-data.rds")
  
  data <- readdata(filename = filename)

  if(is.null(filtermz)){
    mapmoudle(data = data,
              savename = savename,
              mode = mode,
              filtermz = filtermz,
              mapname = mapname,
              ...)
  }else{
    for ( mz in filtermz) {
      pp <- mapmoudle(data = data,
                      savename = savename,
                      mode = mode,
                      filtermz = mz,
                      mapname = paste0(mapname,mz,"-"),
                      ...)
    }
  }
}


#' getpseudotime_plotcelltrajectory
#'
#' 获取拟时序图
#'
#' @param data 数据
#' @param savename 保存名称
#' @param mode 正负离子模式
#' @param savepath 保存路径
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param family 图片字体
#' @param dpi 图片分辨率
#' @param compression tiff格式压缩模式
#' @param text_theme 文本主题
#' @param title_theme 标题主题
#' @param filtermz 绘制mz
#' @param scalemz 是否归一化mz强度
#' @param ... 见[monocle::plot_cell_trajectory()]
#'
#' @export
getpseudotime_plotcelltrajectory  <- function(data,
                                              savename = "all",
                                              mode = "neg",
                                              savepath = "./",
                                              mapname = "",
                                              imagetype = c("png", "pdf"),
                                              other = ggplot2::theme(),
                                              height = 5,
                                              width = 6,
                                              compression = "zip",
                                              dpi = 300,
                                              family = "sans",
                                              filtermz = NULL,
                                              scalemz = T,
                                              pseudotime_plot_moudle = monocle::plot_cell_trajectory,
                                              mapmz = NULL,
                                              seqrange = NULL,
                                              ...) {
  # 加载monocle
  suppressMessages(library("monocle"))

  cds <- data$cds

  if(!is.null(filtermz)){
    intensitydata <- data$spectradata[format(filtermz,nsmall = 5,Trim = T),pData(cds)$Name]
    if(scalemz){
      intensitydata <- scale(intensitydata)
    }
    pData(cds)$intensity <- intensitydata
  }

  if(!is.null(mapmz)){
    cds <- cds[paste0("mz:",format(mapmz,nsmall = 5,Trim = T)),]
  }

  if(!is.null(seqrange)){
    cds <- cds[seqrange,]
  }

  pp <- pseudotime_plot_moudle(cds,...)

  ggplotsave(plot = pp,
             savepath =  savepath,
             family = family,
             mapname = paste0(mapname, savename, "-", mode),
             imagetype = imagetype,
             other = other ,
             height = height,
             width = width,
             compression = compression,
             dpi = dpi)
}


#' getpseudotime_plotcelldensity
#'
#' 获取拟时序密度图
#'
#' @param data 数据
#' @param savename 保存名称
#' @param mode 正负离子模式
#' @param savepath 保存路径
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param family 图片字体
#' @param dpi 图片分辨率
#' @param compression tiff格式压缩模式
#' @param text_theme 文本主题
#' @param title_theme 标题主题
#' @param ... 见[monocle::plot_cell_trajectory()]
#'
#' @export
getpseudotime_plotcelldensity <- function(data,
                                          savename = "all",
                                          mode = "neg",
                                          savepath = "./",
                                          mapname = "",
                                          imagetype = c("png", "pdf"),
                                          other = ggplot2::theme(),
                                          height = 6,
                                          width = 5,
                                          compression = "zip",
                                          dpi = 300,
                                          family = "sans",
                                          filtermz = NULL,
                                          ...) {
  # 加载monocle
  suppressMessages(library("monocle"))

  cds <- data$cds
  df <- pData(cds)

  pp <- ggplot(df,aes(Pseudotime,fill = Cluster,color = Cluster))+
    geom_density(bw = 0.5,alpha = 0.5,size = 0.5)+
    theme_classic()+
    theme(strip.background = element_rect(color = "white", fill = "white"),
          panel.grid = element_blank())

  ggplotsave(plot = pp,
             savepath = savepath,
             family = family,
             mapname = paste0(mapname, savename, "-", mode),
             imagetype = imagetype,
             other = other ,
             height = height,
             width = width,
             compression = compression,
             dpi = dpi)
}

#' getpseudotime_plotcellheatmap
#'
#' 获取拟时序热图
#'
#' @param data 数据
#' @param savename 保存名称
#' @param mode 正负离子模式
#' @param savepath 保存路径
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param family 图片字体
#' @param dpi 图片分辨率
#' @param compression tiff格式压缩模式
#' @param text_theme 文本主题
#' @param title_theme 标题主题
#' @param ... 见[monocle::plot_cell_trajectory()]
#'
#' @export
getpseudotime_plotcellheatmap <- function(data,
                                          savename = "all",
                                          mode = "neg",
                                          savepath = "./",
                                          mapname = "",
                                          imagetype = c("png", "pdf"),
                                          other = ggplot2::theme(),
                                          height = 6,
                                          width = 5,
                                          compression = "zip",
                                          dpi = 300,
                                          family = "sans",
                                          filtermz = NULL,
                                          fillgroup = "Pseudotime",
                                          asp = 1,
                                          ...) {
  # 加载monocle
  suppressMessages(library("monocle"))

  cds <- data$cds
  df <- pData(cds)

  maxx <- max(df$x)+5
  minx <- min(df$x)-5
  maxy <- max(df$y)+5
  miny <- min(df$y)-5

  pp <- ggplot(df,aes(x = x,y = y,fill = get(fillgroup)))+
    geom_tile(width=1,height=1)+
    scale_y_reverse(expand = c(0,0),limits = c(maxy,miny))+
    scale_x_continuous(expand = c(0,0),limits = c(minx,maxx))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.background = element_rect(color = "white", fill = "white"),
          aspect.ratio = (maxy-miny)/(maxx-minx)*asp)+
    xlab("")+
    ylab("")+
    labs(fill = fillgroup)

  ggplotsave(plot = pp,
             savepath = savepath,
             family = family,
             mapname = paste0(mapname, savename, "-", mode),
             imagetype = imagetype,
             other = other ,
             height = height,
             width = width,
             compression = compression,
             dpi = dpi)
}

#' getpseudotime_plotordering
#'
#' 获取拟时序mz筛选图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param color_by 颜色变量名
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_plotordering <- function(mapname = "pseudotime-mzfilter",
                                       height = 5,
                                       width = 5,
                                       pseudotime_plot_moudle = monocle::plot_ordering_genes,
                                       ...){
  getpseudotime_plotcell(mapname = mapname,
                         height = height,
                         width = width,
                         pseudotime_plot_moudle = pseudotime_plot_moudle,
                         ...)
}


#' getpseudotime_Cluster
#'
#' 根据sscc聚类的拟时序图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param color_by 颜色变量名
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_Cluster <- function(mapname = "Pseudotime-Cluster-",
                                  other = list(ggplot2::theme(aspect.ratio = 1),
                                               ggplot2::scale_color_manual(values = SelectColors(palette = "customecol2", n = 50))),
                                  height = 6,
                                  width = 5,
                                  color_by = "Cluster",
                                  ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         color_by = color_by,
                         ...)
}

#' getpseudotime_ClustertoState
#'
#' 根据sscc聚类及State分面的拟时序图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param color_by 颜色变量名
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_ClustertoState <- function(mapname = "Pseudotime-Cluster-StateFacet-",
                                         other = list(ggplot2::theme(aspect.ratio = 1),
                                                      ggplot2::facet_wrap(~State),
                                                      ggplot2::scale_color_manual(values = SelectColors(palette = "customecol2", n = 50))),
                                         height = 15,
                                         width = 10,
                                         color_by = "Cluster",
                                         ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         color_by = color_by,
                         ...)
}

#' getpseudotime_ClustertoSample
#'
#' 根据sscc聚类及Sample分面的拟时序图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param color_by 颜色变量名
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_ClustertoSample <- function(mapname = "Pseudotime-Cluster-SampleFacet-",
                                          other = list(ggplot2::theme(aspect.ratio = 1),
                                                       ggplot2::facet_wrap(~Sample),
                                                       ggplot2::scale_color_manual(values = SelectColors(palette = "customecol2", n = 50))),
                                          height = 15,
                                          width = 10,
                                          color_by = "Cluster",
                                          ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         color_by = color_by,
                         ...)
}

#' getpseudotime_Pseudotime
#'
#' 拟时序图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param color_by 颜色变量名
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_Pseudotime <- function(mapname = "Pseudotime-",
                                     other = list(ggplot2::theme(aspect.ratio = 1),
                                                  ggplot2::scale_colour_gradientn(colours = colorRampPalette(c("white", "blue2", "cyan2", "yellow","brown1", "firebrick3"))(1000))),
                                     height = 6,
                                     width = 5,
                                     color_by = "Pseudotime",
                                     ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         color_by = color_by,
                         ...)
}

#' getpseudotime_PseudotimetoCluster
#'
#' Cluster分面的拟时序图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param color_by 颜色变量名
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_PseudotimetoCluster <- function(mapname = "Pseudotime-ClusterFacet-",
                                              other = list(ggplot2::theme(aspect.ratio = 1),
                                                           ggplot2::facet_wrap(~Cluster),
                                                           ggplot2::scale_colour_gradientn(colours = colorRampPalette(c("white", "blue2", "cyan2", "yellow","brown1", "firebrick3"))(1000))),
                                              height = 15,
                                              width = 10,
                                              color_by = "Pseudotime",
                                              ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         color_by = color_by,
                         ...)
}

#' getpseudotime_PseudotimetoSample
#'
#' Sample分面的拟时序图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param color_by 颜色变量名
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_PseudotimetoSample <- function(mapname = "Pseudotime-SampleFacet-",
                                             other = list(ggplot2::theme(aspect.ratio = 1),
                                                          ggplot2::facet_wrap(~Sample),
                                                          ggplot2::scale_colour_gradientn(colours = colorRampPalette(c("white", "blue2", "cyan2", "yellow","brown1", "firebrick3"))(1000))),
                                             height = 15,
                                             width = 10,
                                             color_by = "Pseudotime",
                                             ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         color_by = color_by,
                         ...)
}

#' getpseudotime_State
#'
#' State的拟时序图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param color_by 颜色变量名
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_State <- function(mapname = "Pseudotime-State-",
                                other = list(ggplot2::theme(aspect.ratio = 1)),
                                height = 6,
                                width = 5,
                                color_by = "State",
                                ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         color_by = color_by,
                         ...)
}

#' getpseudotime_mz
#'
#' mz强度的拟时序图
#'
#' @param filtermz 提供mz
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_mz <- function(filtermz = "105.01784",
                             mapname = "Pseudotime-mz-",
                             other = list(ggplot2::theme(aspect.ratio = 1),
                                          ggplot2::scale_colour_gradientn(colours = colorRampPalette(c("white", "blue2", "cyan2", "yellow","brown1", "firebrick3"))(1000))),
                             height = 6,
                             width = 5,
                             color_by = "intensity",
                             ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         filtermz = filtermz,
                         color_by = color_by,
                         ...)
}

#' getpseudotime_mztoSample
#'
#' mz强度及Sample分面的拟时序图
#'
#' @param filtermz 提供mz
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_mztoSample <- function(filtermz = "105.01784",
                                     mapname = "Pseudotime-mz-SampleFacet-",
                                     other = list(ggplot2::theme(aspect.ratio = 1),
                                                  ggplot2::facet_wrap(~Sample),
                                                  ggplot2::scale_colour_gradientn(colours = colorRampPalette(c("white", "blue2", "cyan2", "yellow","brown1", "firebrick3"))(1000))),
                                     height = 15,
                                     width = 10,
                                     color_by = "intensity",
                                     ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         filtermz = filtermz,
                         color_by = color_by,
                         ...)
}

#' getpseudotime_Clusterdensity
#'
#' 拟时序密度图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param mapmoudle 绘图函数
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_Clusterdensity <- function(mapname = "Pseudotime-Cluster-density-",
                                         other = list(ggplot2::theme(aspect.ratio = 9/16),
                                                      ggplot2::scale_color_manual(values = SelectColors(palette = "customecol2", n = 50)),
                                                      ggplot2::scale_fill_manual(values = SelectColors(palette = "customecol2", n = 50))),
                                         height = 5,
                                         width = 8,
                                         mapmoudle = getpseudotime_plotcelldensity,
                                         ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         mapmoudle = mapmoudle,
                         ...)
}

#' getpseudotime_Clusterdensity
#'
#' Cluster分面的拟时序密度图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param mapmoudle 绘图函数
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_Clusterdensity2 <- function(mapname = "Pseudotime-ClusterFacet-density-",
                                          other = list(ggplot2::theme(aspect.ratio = 9/16),
                                                       ggplot2::facet_wrap(~Cluster),
                                                       ggplot2::scale_color_manual(values = SelectColors(palette = "customecol2", n = 50)),
                                                       ggplot2::scale_fill_manual(values = SelectColors(palette = "customecol2", n = 50))),
                                          height = 10,
                                          width = 15,
                                          mapmoudle = getpseudotime_plotcelldensity,
                                          ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         mapmoudle = mapmoudle,
                         ...)
}

#' getpseudotime_ClusterdensitytoState
#'
#' State分面的拟时序密度图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param mapmoudle 绘图函数
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_ClusterdensitytoState <- function(mapname = "Pseudotime-Cluster-density-StateFacet-",
                                                other = list(ggplot2::theme(aspect.ratio = 9/16),
                                                             ggplot2::facet_wrap(~State),
                                                             ggplot2::scale_color_manual(values = SelectColors(palette = "customecol2", n = 50)),
                                                             ggplot2::scale_fill_manual(values = SelectColors(palette = "customecol2", n = 50))),
                                                height = 10,
                                                width = 15,
                                                mapmoudle = getpseudotime_plotcelldensity,
                                                ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         mapmoudle = mapmoudle,
                         ...)
}

#' getpseudotime_ClusterdensitytoSample
#'
#' Sample分面的拟时序密度图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param mapmoudle 绘图函数
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_ClusterdensitytoSample <- function(mapname = "Pseudotime-Cluster-density-SampleFacet-",
                                                 other = list(ggplot2::theme(aspect.ratio = 9/16),
                                                              ggplot2::facet_wrap(~Sample),
                                                              ggplot2::scale_color_manual(values = SelectColors(palette = "customecol2", n = 50)),
                                                              ggplot2::scale_fill_manual(values = SelectColors(palette = "customecol2", n = 50))),
                                                 height = 10,
                                                 width = 15,
                                                 mapmoudle = getpseudotime_plotcelldensity,
                                                 ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         mapmoudle = mapmoudle,
                         ...)
}

#' getpseudotime_ClusterdensitytoSample2
#'
#' Sample及Cluster分面的拟时序密度图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param mapmoudle 绘图函数
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_ClusterdensitytoSample2 <- function(mapname = "Pseudotime-Cluster&SampleFacet-density-",
                                                  other = list(ggplot2::theme(aspect.ratio = 1),
                                                               ggplot2::facet_grid(vars(Sample),vars(Cluster)),
                                                               ggplot2::scale_color_manual(values = SelectColors(palette = "customecol2", n = 50)),
                                                               ggplot2::scale_fill_manual(values = SelectColors(palette = "customecol2", n = 50))),
                                                  height = 20,
                                                  width = 25,
                                                  mapmoudle = getpseudotime_plotcelldensity,
                                                  ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         mapmoudle = mapmoudle,
                         ...)
}

#' getpseudotime_Clusterexpress
#'
#' 强度分布图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param pseudotime_plot_moudle 拟时序相关绘图函数
#' @param seqrange 绘制mz的范围
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_Clusterexpress <- function(mapname = "Pseudotime-Cluster&express-",
                                         other = list(ggplot2::scale_color_manual(values = SelectColors(palette = "customecol2", n = 50))),
                                         height = 15,
                                         width = 10,
                                         pseudotime_plot_moudle = monocle::plot_genes_in_pseudotime,
                                         color_by = "Cluster",
                                         seqrange = 1:5,
                                         ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         pseudotime_plot_moudle = pseudotime_plot_moudle,
                         seqrange = seqrange,
                         color_by = color_by,
                         ...)
}


#' getpseudotime_Pseudotimeheatmap
#'
#' 拟时序成像图
#'
#' @param mapname 图片名称
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param mapmoudle 绘图函数
#' @param ... 见[getpseudotime_plotcell()]
#'
#' @export
getpseudotime_Pseudotimeheatmap <- function(mapname = "Pseudotime-heatmap-",
                                            other = list(ggplot2::facet_wrap(~Sample),
                                                         ggplot2::scale_fill_gradientn(colours = colorRampPalette(c("white", "blue2", "cyan2", "yellow","brown1", "firebrick3"))(1000))),
                                            height = 10,
                                            width = 10,
                                            mapmoudle = getpseudotime_plotcellheatmap,
                                            ...){
  getpseudotime_plotcell(mapname = mapname,
                         other = other,
                         height = height,
                         width = width,
                         mapmoudle = mapmoudle,
                         ...)
}
