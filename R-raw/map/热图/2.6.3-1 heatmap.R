#!/opt/conda/bin/Rscript

#' 热图可视化
#'
#' @param data 数据
#' @param mapname 保存文件名
#' @param height 图片长度
#' @param width 图片宽度
#' @param log 逻辑，是否log处理
#' @param ... 见[pheatmap()]
#'
#' @export
auto_pheatmap <- function(data,
                          mapname = NA,
                          cluster_rows = T,
                          cluster_cols = F,
                          scale = "auto",
                          treeheight_row = 40 + dim(data)[1] * 0.05,
                          treeheight_col = 40 + dim(data)[2] * 0.05,
                          color = SelectColors(palette = "heatmapcol",n = 255),
                          cellheight = ifelse(600 / dim(data)[1] > 15, 15, 600 / dim(data)[1]),
                          cellwidth = ifelse(600 / dim(data)[2] > 40, 40, 600 / dim(data)[2]),
                          fontsize_row = 0.75 * cellheight,
                          fontsize_col = ifelse(0.5 * cellwidth > 15, 15, 0.5 * cellwidth),
                          angle_col = 90,
                          annotation_row = NA,
                          annotation_col = NA,
                          classfile = "classtype.xlsx",
                          annotation_colors = {if(any(is.na(annotation_col))){NA}else{
                            list(Group = stylefun_group(classfile = classfile,
                                                        object = annotation_col,
                                                        styletype = "fill"))
                          }},
                          border_color = F,
                          show_rownames = T,
                          show_colnames = T,
                          height = (dim(data)[1] * cellheight + length(annotation_col) * 50 + 200 + max(nchar(names(data))) * ifelse(show_colnames, 0.5 * fontsize_col, 0) + ifelse(cluster_cols, treeheight_col, 0)) / 72,
                          width = (dim(data)[2] * cellwidth + length(annotation_row) * 50 + 200 + max(nchar(row.names(data))) * ifelse(show_rownames, 0.5 * fontsize_row, 0) + ifelse(cluster_rows, treeheight_row, 0)) / 72,
                          log = FALSE,
                          savepath = "./",
                          imagetype = c("png","pdf"),
                          dpi = 300,
                          family = "sans",
                          units = "in",
                          ...) {
  # print(annotation_col)
  # print(annotation_colors)
  
  if(any(is.na(data))){
    data[is.na(data)] <- min(as.matrix(data),na.rm = T)
  }
  
  options(warn = -1)
  if(scale=="auto"){
    scale = ifelse(ncol(data)>2,"row","none")
  }
  if (log) {
    if(any(data < 0)){
      warning("数据中有负数，不进行log",immediate. = T)
      data1 <- data
    }else{data1 <- log10(data + 1)}
  }else{data1 <- data}
  
  if (dim(data1)[1] == 1 & cluster_rows) {cluster_rows <- F}
  if (dim(data1)[2] == 1 & cluster_cols) {cluster_cols <- F}
  if (height < 5) {height <- 5}
  if (width < 5) {width <- 5}
  
  plotfile(mapname = mapname,
           height = height,
           width =width,
           savepath = savepath,
           imagetype = imagetype,
           dpi = dpi,
           family = family,
           units = units)
  
  #### heatmap作图函数使用####
  pheatmap::pheatmap(data1,
                     treeheight_row = treeheight_row, # 聚类树的列长度
                     treeheight_col = treeheight_col, # 聚类树行长度
                     scale = scale, # 矩阵有没有进行标准化
                     cluster_cols = cluster_cols, # 按列聚类
                     cluster_rows = cluster_rows, # 按行聚类
                     fontsize_row = fontsize_row, # 行字体大小
                     fontsize_col = fontsize_col, # ；列字体大小
                     cellwidth = cellwidth, cellheight = cellheight, # 设置图片大小
                     color = color, # 设置渐变色
                     border_color = border_color, # 每个小块间是否要用颜色分开
                     annotation_row = annotation_row,
                     annotation_col = annotation_col,
                     height = height,
                     width = width,
                     show_rownames = show_rownames,
                     show_colnames = show_colnames,
                     filename = NA,
                     annotation_colors = annotation_colors,
                     angle_col = angle_col,
                     ...)
  
  plotsave()
  
  return(list(height = height,
              width = width))
}


#' 热图可视化
#'
#' @param data 数据
#' @param name 保存文件名
#' @param type 保存图片格式
#' @param rowgroup 行注释数量
#' @param colgroup 列注释数量
#' @param annotation_row 行注释
#' @param annotation_col 列注释
#' @param ... 见[auto_pheatmap()]
#'
#' @export
auto_heatmap <- function(data,
                         mapname = "heatmap",
                         rowgroup = 0,
                         colgroup = 0,
                         annotation_row = NA,
                         annotation_col = NA,
                         savepath = "./",
                         addgroup = F,
                         ...) {
  
  suppressMessages(library("stringr"))
  #### 数据调取####
  annotation_row1 <- annotation_row
  annotation_col1 <- annotation_col
  
  if (dim(data)[1] < 1+colgroup | dim(data)[2] < 1+rowgroup) {
    
    savetxt(data = "数据为空，未提供热图",
            filename = paste0(savepath,"/说明.txt"))
    
    return()
  }
  
  for(k in 1:nrow(data)){
    row.names(data)[k] <- unlist(strsplit(split = ";\n",x = row.names(data)[k]))[1]
    row.names(data)[k] <- unlist(strsplit(split = "; ",x = row.names(data)[k]))[1]
    if(str_length(row.names(data)[k])>70){
      newname <- paste0(substring(row.names(data)[k],1,35),"...")
      i <- 2
      while (newname %in% row.names(data)) {
        newname <- paste0(substring(row.names(data)[k],1,35),"...-",i)
        i <- i+1
      }
      row.names(data)[k] <- newname
    }
  }
  
  negpos.heatmap <- data
  
  if (rowgroup & colgroup) {
    negpos.heatmap <- data[(1 + colgroup):dim(data)[1], (1 + rowgroup):dim(data)[2], drop = F]
    negpos.heatmap[, ] <- apply(negpos.heatmap, 2, as.numeric)
    annotation_row1 <- data[(1 + colgroup):dim(data)[1], 1:rowgroup, drop = F]
    annotation_col1 <- data[1:colgroup, (1 + rowgroup):dim(data)[2], drop = F]
    annotation_col1 <- data.frame(t(annotation_col1), check.names = F)
    row.names(annotation_col1) <- colnames(negpos.heatmap)
    annotation_row1 <- factorlevel(annotation_row1)
    annotation_col1 <- factorlevel(annotation_col1)
  } else if (rowgroup) {
    negpos.heatmap <- data[, (1 + rowgroup):dim(data)[2], drop = F]
    negpos.heatmap[, ] <- apply(negpos.heatmap, 2, as.numeric)
    annotation_row1 <- data[, 1:rowgroup, drop = F]
    annotation_row1 <- factorlevel(annotation_row1)
  } else if (colgroup) {
    negpos.heatmap <- data[(1 + colgroup):dim(data)[1], , drop = F]
    negpos.heatmap[, ] <- apply(negpos.heatmap, 2, as.numeric)
    annotation_col1 <- data[1:colgroup, , drop = F]
    annotation_col1 <- data.frame(t(annotation_col1), check.names = F)
    row.names(annotation_col1) <- colnames(negpos.heatmap)
    annotation_col1 <- factorlevel(annotation_col1)
  } else {
    negpos.heatmap[, ] <- apply(negpos.heatmap, 2, as.numeric)
  }
  
  returndata <- auto_pheatmap(data = negpos.heatmap,
                              mapname = mapname,
                              annotation_row = annotation_row1,
                              annotation_col = annotation_col1,
                              savepath = savepath,
                              ...)
  
  # print("auto_heatmap运行完成")
  
  return(returndata)
}


#' @export
factorlevel <- function(data) {
  class <- data
  for (i in dim(class)[2]) {
    class[, i] <- factor(x = class[, i], levels = unique(class[, i]))
  }
  return(class)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  map_common_heatmap <- map_autodraw$new(moudle = auto_heatmap,row.names = 1)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "矩阵文件",
                      required = T)
  parser$add_argument("-sh","--sheet",default = NULL,nargs="+",help = "xlsx中的sheet，全部分析请忽略")
  
  # 基本参数
  parser$add_argument("-mn","--mapname", default = NULL, help = "保存文件名")
  parser$add_argument("-s","--savepath",default = "./", help = "保存路径")
  parser$add_argument("-i","--imagetype",default = c("png","pdf"), help = "图片格式",nargs="+",
                      choices = c("jpg","tiff","png","pdf"))
  parser$add_argument("-fa","--family",default = "sans", help = "字体")
  parser$add_argument("-wi","--width",default = 0, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 0, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-cf","--classfile",default = "classtype.xlsx", help = "模板")
  parser$add_argument("-co","--color",default = "heatmapcol",help = "配色设置")
  parser$add_argument("-lg","--log", default = F, action = "store_true",help="是否log处理")
  parser$add_argument("-cr","--cluster_rows", default = F, action = "store_true", help = "按行聚类")
  parser$add_argument("-cc","--cluster_cols", default = F, action = "store_true", help = "按列聚类")
  parser$add_argument("-sl","--scale", default = "auto", help = "数据是否（按行、列）进行标准化，默认自动判断",
                      choices = c("auto","row","column","none"))
  parser$add_argument("-bc","--border_color", default = F, action = "store_true", help = "每个格子是否要用颜色分开")
  parser$add_argument("-sr","--show_rownames", default= F, action = "store_true", help = "显示行名")
  parser$add_argument("-sc","--show_colnames", default= F, action = "store_true", help = "显示列名")
  parser$add_argument("-rg","--rowgroup", default = 0,type = "integer",help = "行注释数量")
  parser$add_argument("-cg","--colgroup", default = 0,type = "integer",help = "列注释数量")
  # parser$add_argument("-ar","--annotation_row", default = F, help = "行注释")
  # parser$add_argument("-ac","--annotation_col", default = F, help = "列注释")
  
  args <- parser$parse_args()
  
  args$color <- SelectColors(palette = args$color,n = 255)
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  heatmapargs <- do.call(what = map_common_heatmap, args = args)
}

#' 根据文件进行heatmap可视化
#' 
#' @export
map_common_heatmap <- map_autodraw$new(moudle = auto_heatmap,row.names = 1)$draw
