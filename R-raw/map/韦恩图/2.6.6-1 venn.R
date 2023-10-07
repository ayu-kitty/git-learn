#!/opt/conda/bin/Rscript

#' 绘制韦恩图、upset、花瓣图
#'
#' @param data 数据
#' @param name 保存文件名
#' @param type 保存图片格式
#' @param height 图片长度
#' @param width 图片宽度
#' @param ppi 图片分辨率
#' @param family 图片字体
#' @param savedata 逻辑，是否保存数据
#' @param xlsxname 保存venn结果文件名
#' @param upset 逻辑，是否使用upset绘图
#' @param ... 见[upset()]或[venn.diagram()]或[flower_plot()]
#'
#' @export
auto_venn <- function(data,
                      mapname = "Venn",
                      width = 12,
                      height = 12,
                      savepath = "./",
                      imagetype = c("png","pdf"),
                      dpi = 300,
                      family = "sans",
                      units = "in",
                      cat.col = "black",
                      col = rep("transparent",length(lengths(data))),
                      fillpalette = "venncol2",
                      fill = SelectColors(palette = fillpalette,n = length(lengths(data))),
                      alpha = 0.40,
                      cex = 1+1/length(lengths(data)),
                      cat.cex = 1 + 10 / (max(nchar(names(data)))*length(lengths(data))),
                      cat.fontface = "plain",
                      fontface = "plain",
                      margin = 0.25,
                      scale = F,
                      ext.text = F,
                      euler.d = F,
                      savedata = T,
                      xlsxname = paste0(mapname,".xlsx"),
                      upset = F,
                      nsets = length(data),
                      nintersects = NA,
                      text.scale = c(1.3, 1.3, 1, 1, 1.3, 1),
                      order.by = "freq",
                      ...) {
  suppressWarnings(library(scales))
  
  rawdata <- data
  
  data1 <- Reduce(union, rawdata)
  data3 <- data.frame(ID = data1)
  i <- 1
  while (i <= length(data)) {
    # unlist(data[i])
    data2 <- data.frame(c(rep(0, length(data1))))
    data2[data3[, 1] %in% unlist(data[i]), ] <- 1
    names(data2) <- names(data[i])
    data3 <- data.frame(data3, data2, check.names = F)
    i <- i + 1
  }
  
  if (savedata) {
    savexlsx1(data = data3, 
              filename= paste0(savepath,"/",xlsxname),
              sheet = "venn")
  }
  
  if (upset | length(data) > 5) {
    
    if(upset){
      num <- ncol(data3)
      
      if(num < 6){
        width <- 8
        height <- 6
      }else if(num < 11){
        width <- 16
        height <- 8
      }else if(num < 20){
        width <- 32
        height <- 10
      }else{
        width <- 48
        height <- 12
      }
      
      plotfile(mapname = mapname,
               height = height,
               width =width,
               savepath = savepath,
               imagetype = imagetype,
               dpi = dpi,
               family = family,
               units = units)
      
      print(UpSetR::upset(data = data3,
                          nsets = nsets,
                          nintersects = nintersects,
                          text.scale = text.scale,
                          order.by = order.by,
                          sets.bar.color = alpha(colour = fill,0.4),
                          ...))
    }else{
      plotfile(mapname = mapname,
               height = height,
               width =width,
               savepath = savepath,
               imagetype = imagetype,
               dpi = dpi,
               family = family,
               units = units)
      
      flower_plot(data = rawdata, ...)
    }

    plotsave()
  } else {
    num <- length(data)
    if(num == 5){
      cat.pos <- c(0,335,220,180,50)
      cat.dist <- rep(0.2, 5)
    }else if(num == 4){
      cat.pos <- c(-30, 30, 0, 0)
      cat.dist <- c(0.3,0.3, 0.15, 0.15)
    }else if(num == 3){
      cat.pos <- c(-30, 30, 180)
      cat.dist <- c(0.1, 0.1, 0.05)
    }else if(num == 2){
      cat.pos <- c(-30, 30)
      cat.dist <- rep(0.1, 2)
    }else if(num == 1){
      cat.pos <- c(0)
      cat.dist <- 0.025
    }
    
    venn <- VennDiagram::venn.diagram(rawdata,
                                      resolution = resolution,
                                      filename = NULL,
                                      col = col,
                                      fill = fill,
                                      alpha = alpha,
                                      cex = cex,
                                      cat.col = cat.col,
                                      cat.cex = cat.cex,
                                      cat.dist = cat.dist,
                                      cat.fontface = cat.fontface,
                                      margin = margin,
                                      scale = scale,
                                      ext.text = ext.text,
                                      euler.d = euler.d,
                                      main.fontfamily = family,
                                      sub.fontfamily = family,
                                      fontfamily = family,
                                      cat.fontfamily = family,
                                      cat.pos = cat.pos,
                                      fontface = fontface,
                                      ...)
    unlink(list.files(pattern = "^VennDiagram.*\\.log"))
    ggplotsave(plot = venn,
               mapname = mapname,
               height = height,
               width =width,
               savepath = savepath,
               imagetype = imagetype,
               dpi = dpi,
               family = family,
               units = units)
  }
  
}

#' @export
flower_plot <- function(data = data,
                        ellipse_col = NULL,
                        circle_col = rgb(255, 0, 0, 240, maxColorValue = 255)) {
  
  sample <- names(data)
  value <- lengths(data)
  mid <- length(Reduce(intersect, data))
  
  par(bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1, 1, 1, 1))
  plot(c(0, 15), c(0, 15), type = "n")
  n <- length(sample)
  deg <- 360 / n
  
  if (n %% 6 == 0) {
    n1 <- as.vector(matrix(data = c(1:n), nrow = 6, byrow = T, dimnames = NULL))
  } else if (n < 12) {
    n2 <- ceiling(n / 6)
    n3 <- n %% 6
    n4 <- matrix(data = c(1:(n2 * n3)), ncol = n2, byrow = T, dimnames = NULL)
    n5 <- matrix(data = c((n2 * n3 + 1):n), ncol = n2 - 1, byrow = T, dimnames = NULL)
    n1 <- c(n4[, -n2], n5, n4[, n2])
  } else {
    n2 <- ceiling(n / 6)
    n3 <- n %% 6
    n4 <- matrix(data = c(1:(n2 * n3)), ncol = n2, byrow = T, dimnames = NULL)
    n5 <- matrix(data = c((n2 * n3 + 1):n), ncol = n2 - 1, byrow = T, dimnames = NULL)
    n1 <- c(as.vector(rbind(n4[, -n2], n5)), n4[, n2])
  }
  
  t <- 1
  while (t <= n) {
    if (t %% 6 == 1) {
      col <- rgb(200, (t %/% 6) * 200 / (n %/% 6 + 1), 0, 100, maxColorValue = 255)
    } else if (t %% 6 == 2) {
      col <- rgb(200 - (t %/% 6) * 200 / (n %/% 6 + 1), 200, 0, 100, maxColorValue = 255)
    } else if (t %% 6 == 3) {
      col <- rgb(0, 200, (t %/% 6) * 200 / (n %/% 6 + 1), 100, maxColorValue = 255)
    } else if (t %% 6 == 4) {
      col <- rgb(0, 200 - (t %/% 6) * 200 / (n %/% 6 + 1), 200, 100, maxColorValue = 255)
    } else if (t %% 6 == 5) {
      col <- rgb((t %/% 6) * 200 / (n %/% 6 + 1), 0, 200, 100, maxColorValue = 255)
    } else if (t %% 6 == 0) {
      col <- rgb(200, 0, 200 - (t %/% 6 - 1) * 200 / (n %/% 6), 100, maxColorValue = 255)
    }
    
    if (!is.null(ellipse_col)) {
      col <- ellipse_col[n1[t] - 1]
    }
    
    plotrix::draw.ellipse(x = 7.5 + 3 * 0.8 * (1 - 2 / n) * sin(deg * (n1[t] - 1) * pi / 180),
                          y = 7.5 + 3 * 0.8 * (1 - 1.5 / n) * cos(deg * (n1[t] - 1) * pi / 180),
                          col = col,
                          border = col,
                          a = 0.2 + 6 / n, b = 2.5, angle = -deg * (n1[t] - 1))
    text(x = 7.5 + 4.8 * (1 - 2 / n) * sin(deg * (n1[t] - 1) * pi / 180),
         y = 7.5 + 4.8 * (1 - 2 / n) * cos(deg * (n1[t] - 1) * pi / 180),
         value[n1[t]], cex = 0.4 + 5 / n, font = 2)
    
    if (n >= 13) {
      if (deg * (n1[t] - 1) <= 180 && deg * (n1[t] - 1) >= 0) {
        srt <- -deg * (n1[t] - 1) + 90
        adj <- 0
      } else {
        srt <- -deg * (n1[t] - 1) - 90
        adj <- 1
      }
    } else {
      if (deg * (n1[t] - 1) < 180 && deg * (n1[t] - 1) > 0) {
        srt <- 0
        adj <- 0
      } else if (deg * (n1[t] - 1) == 180 || deg * (n1[t] - 1) == 0) {
        srt <- 0
        adj <- 0.5
      } else {
        srt <- 0
        adj <- 1
      }
    }
    
    text(x = 7.5 + 5.6 * (1 - 1.2 / n) * sin(deg * (n1[t] - 1) * pi / 180),
         y = 7.5 + 5.6 * (1 - 1.2 / n) * cos(deg * (n1[t] - 1) * pi / 180),
         sample[n1[t]],
         srt = srt, adj = adj, cex = 0.5 + 3 / n, font = 2)
    
    t <- t + 1
  }
  plotrix::draw.circle(x = 7.5, y = 7.5, r = 0.14 + 5 / n, col = circle_col, border = circle_col)
  text(x = 7.5, y = 7.5, mid, adj = 0.5, cex = 0.5 + 6 / n, font = 2)
}


#' 绘制韦恩图、upset、花瓣图
#'
#' @param data 数据
#' @param name 保存文件名
#' @param type 保存图片格式
#' @param ... 见[auto_venn()]
#'
#' @export
auto_vennmap <- function(data,
                         mapname = "Venn",
                         selectcol = ifelse(is.list(data),colnames(data[[1]])[1],""),
                         savepath = "./",
                         savediff_filter = T,
                         savediff_filtername = paste0(mapname,".xlsx"),
                         savedata = T,
                         ...) {
  # data2 <<- data
  if (is.data.frame(data)) {
    rawdata <- as.list(data)
  } else {
    if ("diff_filter" %in% unlist(lapply(data,class))){
      rawdata <- NULL
      infodata <- data[[1]]$args$info
      infodata <- infodata[row.names(data[[1]]$result$diffdata),]
      infodata2 <- infodata[,0,drop = F]
      diffdata <- infodata[,0,drop = F]
      expressdata <- infodata[,0,drop = F]
      
      for (i in 1:length(data)) {
        rawdata1 <- list(num = row.names(data[[i]]$result$filterdata))
        groupname <- gsub(pattern = "/",replacement = "-vs-",x = data[[i]]$args$compare)
        names(rawdata1) <- groupname 
        rawdata <- c(rawdata,rawdata1)
        infodata2[,groupname ] <- 0
        infodata2[row.names(data[[i]]$result$filterdata),groupname ] <- 1
        expressdata2 <- data[[i]]$args$data
        expressdata2 <- expressdata2[,!(colnames(expressdata2) %in% colnames(expressdata)),drop = F]
        expressdata <- merge(expressdata,expressdata2,by = 0)
        row.names(expressdata) <- expressdata[,1]
        expressdata <- expressdata[,-1,drop = F]
        
        diffdata2 <- data[[i]]$result$diffdata
        diffdata2 <- diffdata2[,!(colnames(diffdata2) %in% colnames(diffdata))]
        
        if("FoldChange" %in% colnames(diffdata2)){
          colnames(diffdata2)[colnames(diffdata2) == "FoldChange"] <- paste0("FoldChange(",groupname,")")
        }
        if("log2FoldChange" %in% colnames(diffdata2)){
          colnames(diffdata2)[colnames(diffdata2) == "log2FoldChange"] <- paste0("log2FoldChange(",groupname,")")
        }
        if("VIP" %in% colnames(diffdata2)){
          colnames(diffdata2)[colnames(diffdata2) == "VIP"] <- paste0("VIP(",groupname,")")
        }
        if("Regulation" %in% colnames(diffdata2)){
          # diffdata2 <- diffdata2[,colnames(diffdata2) != "Regulation"]
          colnames(diffdata2)[colnames(diffdata2) == "Regulation"] <- paste0("Regulation(",groupname,")")
        }
        if("p-value" %in% colnames(diffdata2)){
          colnames(diffdata2)[colnames(diffdata2) == "p-value"] <- paste0("p-value(",groupname,")")
        }
        if("q-value" %in% colnames(diffdata2)){
          colnames(diffdata2)[colnames(diffdata2) == "q-value"] <- paste0("q-value(",groupname,")")
        }
        
        diffdata <- merge(diffdata,diffdata2,by = 0,sort = F)
        row.names(diffdata) <- diffdata[,1]
        diffdata <- diffdata[,-1,drop = F]
      }
      diffdata <- diffdata[,order(colnames(diffdata)),drop = F]
      
      infodata3 <- cbind(infodata[,1,drop = F],infodata2)
      infodata <- cbind(infodata3,infodata[,2:dim(infodata)[2],drop = F])
      infodata <- merge(infodata,diffdata,by = 0,sort = F)
      row.names(infodata) <- infodata[,1]
      infodata <- infodata[,-1,drop = F]
      
      expressdata <- expressdata[,!duplicated(colnames(expressdata)),drop = F]
      expressdata <- expressdata[,order(colnames(expressdata)),drop = F]
      infodata <- merge(infodata,expressdata,by = 0,sort = F)
      row.names(infodata) <- infodata[,1]
      infodata <- infodata[,-1,drop = F]
      
      intersectiondata <- infodata[row.names(infodata2)[apply(infodata2, 1, function(x){all(x == 1)})],]
      uniondata <- infodata[row.names(infodata2)[apply(infodata2, 1, function(x){any(x == 1)})],]
      
      if(savediff_filter){
        wb <- openxlsx::createWorkbook()
        wb <- addsheet(data = infodata,wb = wb,sheet = "venn")
        wb <- addsheet(data = intersectiondata,wb = wb,sheet = "交集")
        wb <- addsheet(data = uniondata,wb = wb,sheet = "并集")
        savewb(wb = wb,
               filename = paste0(savepath,"/",savediff_filtername),
               overwrite = T)
        savedata <- F
      }
      
    }else{
      rawdata <- data
    }
  }
  
  for (i in 1:length(data)) {
    if(is.data.frame(data[[i]])){
      rawdata[[i]] <- data[[i]][, selectcol]
    }
    rawdata[[i]] <- rawdata[[i]][!is.na(rawdata[[i]])]
  }
  
  rawdata <- rawdata[order(lengths(rawdata),decreasing = T)]
  # rawdata <<- rawdata
  
  auto_venn(data = rawdata,
            savepath = savepath,
            savedata = savedata,
            mapname = mapname,
            ...)
  
  return("auto_vennmap运行完成")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  options(warn = -1)
  
  suppressMessages(library("lmbio"))
  
  map_common_venn <- map_autodraw$new(moudle = auto_vennmap)$draw
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",nargs="+",
                      help = "火山图矩阵文件",required = T)
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
  parser$add_argument("-u","--upset", default = F, action = "store_true", help="是否绘制upset")
  
  args <- parser$parse_args()
  
  if(args$width == 0){ args$width <- NULL}
  if(args$height == 0){ args$height <- NULL}
  
  if(length(args$filename) > 1){
    args$filename <- as.list(args$filename)
  }else{
    args$sheet <- list(args$sheet)
  }
  
  # print(args)
  
  result <- do.call(what = map_common_venn,args = args) 
  
}

#' 根据文件进行venn可视化
#' 
#' @export
map_common_venn <- map_autodraw$new(moudle = auto_vennmap)$draw
