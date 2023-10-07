#!/opt/conda/bin/Rscript

#' getdimredmap
#'
#' 获取umap和tsne相关的图
#'
#' @param datapath sscc数据路径
#' @param samplename 样本名称
#' @param savename 保存名称
#' @param mode 正负离子模式
#' @param savepath 保存路径
#' @param moudle 可视化模块
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param ssccmap 逻辑，是否根据sscc绘制
#' @param ssccpath sscc数据路径
#' @param ssccname sscc名称
#' @param intensitymap 逻辑，是否根据强度绘制
#' @param mapmz 绘制的mz
#' @param annotation_colors 颜色
#' @param annotation_fill 颜色
#' @param annotation_shape 形状
#' @param annosamplename 逻辑，是否根据样本绘制
#' @param areamap 逻辑，是否根据选区绘制
#' @param areaname 选区名
#' @param areapath 选区数据路径
#' @param noareaname 无区域名称
#' @param noareacolor 无区域颜色
#' @param saveareaname 保存区域名
#' @param ... 见moudle函数
#'
#' @export
getdimredmap <- function(datapath = "./sample/cluster/tsne/",
                         samplename = NULL,
                         savename = samplename,
                         mode = "neg",
                         savepath = "./",
                         ssccmap = T,
                         ssccpath = "./sample/cluster/sscc/",
                         ssccname = samplename,
                         intensitymap = F,
                         mapmz = 1,
                         moudle = NULL,
                         mapname = "",
                         imagetype = c("png", "pdf"),
                         qualitativepath = NULL,
                         annotation_colors = list(Cluster = SelectColors(palette = "customecol2",n=50),
                                                  Sample = SelectColors(palette = "cold",n = 30),
                                                  Area = SelectColors(palette = "blindless",n = 30),
                                                  Group = SelectColors(palette = "blindless",n = 30)),
                         annotation_fill = list(Cluster = SelectColors(palette = "customecol2",n=50),
                                                Sample = SelectColors(palette = "cold",n = 30),
                                                Area = SelectColors(palette = "blindless",n = 30),
                                                Group = SelectColors(palette = "blindless",n = 30)),
                         annotation_shape = rep(21,50),
                         annosamplename = T,
                         areamap = F,
                         areaname = "ALL",
                         areapath = "./area/data/",
                         noareaname = "none",
                         noareacolor = "grey",
                         saveareaname = areaname[1],
                         groupmap = F,
                         groupname = "ALL",
                         grouppath = "项目登记单.xlsx",
                         groupdatapath = areapath,
                         nogroupname = "none",
                         nogroupcolor= "grey",
                         savegroupname = groupname[1],
                         stoprdsnoexist = T,
                         filename = paste0(datapath,"/",samplename,"-", mode, "-data.rds"),
                         ssccrds = paste0(ssccpath,ssccname,"-", mode, "-data.rds"),
                         ...){
  suppressMessages(library("Cardinal"))
  print(datapath)
  # rds <- paste0(datapath,samplename,"-", mode, "-data.rds")
  data <- readdata(filename = filename)
  
  plotdata <- data$plotdata
  getmapname <- colnames(plotdata)[1:2]
  
  if(ssccmap){
    # ssccrds <- paste0(ssccpath,ssccname,"-", mode, "-data.rds")
    ssccdata <- readdata(filename = ssccrds)
    ssccplotdata <- ssccdata$plotdata
    
    plotdata <- merge(plotdata,ssccplotdata[,c("Name","Cluster")],by = "Name")
    
    if(is.function(moudle)){
      plotdata2 <- plotdata[,c("Name",getmapname,"Cluster","Sample")]
      plotdata2 <- plotdata2[order(plotdata2$Cluster), ]
      num <- length(unique(plotdata2$Cluster))
      fill <- annotation_fill$Cluster[1:num]
      colour <- annotation_colors$Cluster[1:num]
      shape <- annotation_shape[1:num]
      plotdata2[,"Group"] <- plotdata2[,"Cluster"]
      returndata <- moudle(data = plotdata2,
                           mapname = paste0(mapname,savename,"-", mode),
                           savepath = savepath,
                           imagetype = imagetype,
                           point_fill = fill,
                           point_colour = colour,
                           point_shape = shape,
                           ...)
    }
    
    if(annosamplename & length(unique(plotdata$Sample)) > 1){
      plotdata2 <- plotdata[,c("Name",getmapname,"Sample")]
      plotdata2 <- plotdata2[order(plotdata2$Sample), ]
      num <- length(unique(plotdata2$Sample))
      fill <- annotation_fill$Sample[1:num]
      colour <- annotation_colors$Sample[1:num]
      shape <- annotation_shape[1:num]
      plotdata2[,"Group"] <- plotdata2[,"Sample"]
      returndata <- moudle(data = plotdata2,
                           mapname = paste0(mapname,savename,"-sample-", mode),
                           savepath = savepath,
                           imagetype = imagetype,
                           point_fill = fill,
                           point_colour = colour,
                           point_shape = shape,
                           ...)
      
      plotdata2 <- plotdata[,c("Name",getmapname,"Cluster","Sample")]
      plotdata2 <- plotdata2[order(plotdata2$Cluster), ]
      num <- length(unique(plotdata2$Cluster))
      fill <- annotation_fill$Cluster[1:num]
      colour <- annotation_colors$Cluster[1:num]
      shape <- annotation_shape[1:num]
      plotdata2[,"Group"] <- plotdata2[,"Cluster"]
      returndata <- moudle(data = plotdata2,
                           mapname = paste0(mapname,savename,"-sample-facet-", mode),
                           savepath = savepath,
                           imagetype = imagetype,
                           addfacet = T,
                           point_fill = fill,
                           point_colour = colour,
                           point_shape = shape,
                           ...)
    }
  }
  
  if(intensitymap){
    if(is.null(mapmz)){
      mzdata <- data$spectradata
    }else{
      mzdata <- data$spectradata[mapmz,]
    }
    
    if(dim(mzdata)[1]>0){
      mzname <- rownames(mzdata)
      mzdata <- as.data.frame(t(mzdata))
      plotdata <- merge(plotdata,mzdata,by.x = "Name",by.y=0)
      
      for ( i in 1:length(mzname)) {
        plotdata2 <- plotdata[,c("Name",getmapname,mzname[i],"Sample")]
        plotdata2[,"Group"] <- plotdata2[,mzname[i]]
        returndata <- moudle(data = plotdata2,
                             mapname = paste0(mapname,savename,"-",mzname[i],"-", mode),
                             savepath = savepath,
                             imagetype = imagetype,
                             gradientmode = T,
                             title = mzname[i],
                             ...)
        
        if(annosamplename & length(unique(plotdata$Sample)) > 1){
          returndata <- moudle(data = plotdata2,
                               mapname = paste0(mapname,savename,"-",mzname[i],"-facet-", mode),
                               savepath = savepath,
                               imagetype = imagetype,
                               gradientmode = T,
                               title = mzname[i],
                               addfacet = T,
                               ...)
        }
      }
    }
  }
  
  if(areamap){
    plotdata[,"Area"] <- noareaname
    
    for ( i in 1:length(areaname)) {
      areafile <- list.files(path = areapath,
                             pattern = paste0("-",areaname[i],"-",mode,"-name.rds$"),
                             full.names = T)
      for ( j in 1:length(areafile)) {
        print(paste0(areafile[j],"读取中"))
        if(file.exists(areafile[j])){
          areadata <- readRDS(areafile[j])
        }else{
          if(stoprdsnoexist){
            stop(paste0(areafile[j],"不存在"))
          }else{
            warning(paste0(areafile[j],"不存在"),immediate. = T)
            next
          }
        }
        plotdata[plotdata$Name %in% areadata,"Area"] <- areaname[i]
      }
    }
    
    if(is.function(moudle)){
      plotdata2 <- plotdata[,c("Name",getmapname,"Area","Sample")]
      plotdata2$Area <- factor(plotdata2$Area,levels=c(noareaname,areaname))
      plotdata2 <- plotdata2[order(plotdata2$Area),]
      
      num <- length(unique(plotdata2$Area))
      fill <- c(if(any(plotdata2$Area==noareaname)){noareacolor}else{NULL},annotation_fill$Area)[1:num]
      colour <- c(if(any(plotdata2$Area==noareaname)){noareacolor}else{NULL},annotation_colors$Area)[1:num]
      shape <- annotation_shape[1:num]
      
      plotdata2[,"Group"] <- plotdata2[,"Area"]
      returndata <- moudle(data = plotdata2,
                           mapname = paste0(mapname,savename,"-area-",saveareaname,"-", mode),
                           savepath = savepath,
                           type = imagetype,
                           point_fill = fill,
                           point_colour = colour,
                           point_shape = shape,
                           ...)
      
      if(annosamplename & length(unique(plotdata$Sample)) > 1){
        returndata <- moudle(data = plotdata2,
                             mapname = paste0(mapname,savename,"-area-",saveareaname,"-facet-", mode),
                             savpath = savepath,
                             imagetype = imagetype,
                             point_fill = fill,
                             point_colour = colour,
                             point_shape = shape,
                             addfacet = T,
                             ...)
      }
    }
  }
  
  if(groupmap){
    plotdata[,"Group"] <- nogroupname
    groupdata <- readdata(filename = grouppath,sheet = "分组信息")
    
    for ( i in 1:length(groupname)) {
      
      getgroupdata <- groupdata[groupdata$分组 == groupname[i] & (groupdata$模式 == mode | groupdata$模式 =="both"),]
      
      areafile <- paste0(groupdatapath,getgroupdata$样品,"-",getgroupdata$选区,"-",mode,"-name.rds")
      
      for ( j in 1:length(areafile)) {
        print(paste0(areafile[j],"读取中"))
        if(file.exists(areafile[j])){
          areadata <- readRDS(areafile[j])
        }else{
          if(stoprdsnoexist){
            stop(paste0(areafile[j],"不存在"))
          }else{
            warning(paste0(areafile[j],"不存在"),immediate. = T)
            next
          }
        }
        
        plotdata[plotdata$Name %in% areadata,"Group"] <- groupname[i]
      }
    }
    
    if(is.function(moudle)){
      plotdata2 <- plotdata[,c(getmapname,"Group","Sample")]
      plotdata2$Group <- factor(plotdata2$Group,levels=c(nogroupname,groupname))
      plotdata2 <- plotdata2[order(plotdata2$Group),]
      
      num <- length(unique(plotdata2$Group))
      fill <- c(if(any(plotdata2$Group==nogroupname)){nogroupcolor}else{NULL},annotation_fill$Group)[1:num]
      colour <- c(if(any(plotdata2$Group==nogroupname)){nogroupcolor}else{NULL},annotation_colors$Group)[1:num]
      shape <- annotation_shape[1:num]
      
      returndata <- auto_draw(data = plotdata2,
                              name = paste0(mapname,savename,"-group-",savegroupname,"-", mode),
                              moudle = moudle,
                              frontname = savepath,
                              type = imagetype,
                              fill = fill,
                              colour = colour,
                              shape = shape,
                              ...)
      
      if(annosamplename & length(unique(plotdata$Sample)) > 1){
        returndata <- auto_draw(data = plotdata2,
                                name = paste0(mapname,savename,"-group-",savegroupname,"-facet-", mode),
                                moudle = moudle,
                                frontname = savepath,
                                type = imagetype,
                                fill = fill,
                                colour = colour,
                                shape = shape,
                                addfacet = T,
                                ...)
      }
    }
  }
  
  savetxt(data = plotdata,
          filename = paste0(savepath,"/",mapname,mode, ".csv"),
          row.names = F,
          sep = ",")
}

#' gettsnemap
#'
#' 绘制tsne相关图
#'
#' @param height 图片长度
#' @param width 图片宽度
#' @param datapath 数据路径
#' @param moudle 绘图函数
#' @param mapname 图片名称
#' @param stat_ellipse 绘制置信区间
#' @param y.axis y轴0位置添加线
#' @param x.axis x轴0位置添加线
#' @param point_size 点大小
#' @param annotation_colors 颜色
#' @param annotation_fill 颜色
#' @param annotation_shape 形状
#' @param ... 见[getdimredmap()]
#'
#' @export
gettsnemap <- function(datapath = "./sample/cluster/tsne/",
                       moudle = mulstatisticsmap,
                       mapname = "Tsne-",
                       stat_ellipse = ggplot2::theme(),
                       hline = ggplot2::theme(),
                       vline = ggplot2::theme(),
                       point_size = 1,
                       width = 8,
                       height = 7,
                       annotation_colors = list(Cluster = SelectColors(palette = "customecol2",n=50),
                                                Sample = SelectColors(palette = "cold",n = 30),
                                                Area = SelectColors(palette = "blindless",n = 30),
                                                Group = SelectColors(palette = "blindless",n = 30)),
                       annotation_fill = list(Cluster = SelectColors(palette = "customecol2",n=50),
                                              Sample = SelectColors(palette = "cold",n = 30),
                                              Area = SelectColors(palette = "blindless",n = 30),
                                              Group = SelectColors(palette = "blindless",n = 30)),
                       annotation_shape = rep(21,50),
                       ...){
  getdimredmap(datapath = datapath,
               moudle = moudle,
               mapname = mapname,
               stat_ellipse = stat_ellipse,
               point_size = point_size,
               width = width,
               height = height,
               annotation_colors = annotation_colors,
               annotation_fill = annotation_fill,
               annotation_shape = annotation_shape,
               hline = hline,
               vline = vline,
               ...)
}

#' getumapmap
#'
#' 绘制umap相关图
#'
#' @param height 图片长度
#' @param width 图片宽度
#' @param datapath 数据路径
#' @param moudle 绘图函数
#' @param mapname 图片名称
#' @param stat_ellipse 绘制置信区间
#' @param y.axis y轴0位置添加线
#' @param x.axis x轴0位置添加线
#' @param point_size 点大小
#' @param annotation_colors 颜色
#' @param annotation_fill 颜色
#' @param annotation_shape 形状
#' @param ... 见[getdimredmap()]
#'
#' @export
getumapmap <- function(datapath = "./sample/cluster/umap/",
                       moudle = mulstatisticsmap,
                       mapname = "Umap-",
                       stat_ellipse = ggplot2::theme(),
                       hline = ggplot2::theme(),
                       vline = ggplot2::theme(),
                       point_size = 1,
                       width = 8,
                       height = 7,
                       annotation_colors = list(Cluster = SelectColors(palette = "customecol2",n=50),
                                                Sample = SelectColors(palette = "cold",n = 30),
                                                Area = SelectColors(palette = "blindless",n = 30),
                                                Group = SelectColors(palette = "blindless",n = 30)),
                       annotation_fill = list(Cluster = SelectColors(palette = "customecol2",n=50),
                                              Sample = SelectColors(palette = "cold",n = 30),
                                              Area = SelectColors(palette = "blindless",n = 30),
                                              Group = SelectColors(palette = "blindless",n = 30)),
                       annotation_shape = rep(21,50),
                       ...){
  getdimredmap(datapath = datapath,
               moudle = moudle,
               mapname = mapname,
               stat_ellipse = stat_ellipse,
               point_size = point_size,
               width = width,
               height = height,
               annotation_colors = annotation_colors,
               annotation_fill = annotation_fill,
               annotation_shape = annotation_shape,
               hline = hline,
               vline = vline,
               ...)
}
