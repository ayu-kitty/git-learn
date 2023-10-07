#' 选区质谱成像图
#'
#' @param filename 文件路径
#' @param areafile 区域文件路径
#' @param savepath 保存路径
#' @param savename 保存名称
#' @param type 图片格式
#' @param width 图片宽度
#' @param height 图片长度
#' @param asp 图像长宽比
#' @param lightmode 成像模式
#' @param fun 处理函数
#' @param smooth.image 是否平滑处理
#' @param superpose 是否分面绘制
#' @param normalize.image 归一化模式
#' @param contrast.enhance 成像强对比
#' @param xlab x轴标签
#' @param ylab y轴标签
#' @param color 成像颜色
#' @param family 字体
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param ... 见[Cardinal::image()]
#'
#' @export
imzmlareaimage <- function(filename,
                           areafile,
                           savepath = "./",
                           savename = "test.png",
                           type = "png",
                           width = 8,
                           height = 7,
                           asp = 1,
                           lightmode = F,
                           fun = mean,
                           smooth.image = "gaussian",
                           superpose = F,
                           normalize.image = "linear",
                           contrast.enhance = "none",
                           xlab = "",
                           ylab = "",
                           color = c("blue2", "cyan2", "yellow",
                                     "brown1", "firebrick3"),
                           family = "sans",
                           mass.range = NULL,
                           resolution = 5,
                           units = "ppm",
                           ...) {
  print("质谱成像绘图开始")
  
  # 判断path是否存在，不存在进行创建
  if (!dir.exists(savepath)) {
    dir.create(path = savepath, recursive = T)
  }
  
  suppressMessages(library("Cardinal"))
  
  mse <- readMSIData(filename,
                     attach.only = TRUE,
                     mass.range = mass.range, resolution = resolution, units = units)
  
  samplearea <- readRDS(file = areafile)
  
  tic <- pixelApply(mse, fun)
  
  xmax <- max(coord(mse)$x)
  xmin <- min(coord(mse)$x)
  ymax <- max(coord(mse)$y)
  ymin <- min(coord(mse)$y)
  xlength <- xmax - xmin + 1
  ylength <- ymax - ymin + 1
  
  areaframe <- data.frame(
    x = coord(mse)$x,
    y = coord(mse)$y,
    x_y = paste0(coord(mse)$x, "_", coord(mse)$y),
    area = samplearea,
    frame = F,
    up = paste0(coord(mse)$x, "_", coord(mse)$y + 1),
    down = paste0(coord(mse)$x, "_", coord(mse)$y - 1),
    left = paste0(coord(mse)$x - 1, "_", coord(mse)$y),
    right = paste0(coord(mse)$x + 1, "_", coord(mse)$y),
    upleft = paste0(coord(mse)$x - 1, "_", coord(mse)$y + 1),
    upright = paste0(coord(mse)$x + 1, "_", coord(mse)$y + 1),
    downleft = paste0(coord(mse)$x - 1, "_", coord(mse)$y - 1),
    downright = paste0(coord(mse)$x + 1, "_", coord(mse)$y - 1)
  )
  
  for (i in 1:dim(areaframe)[1]) {
    if (!all(areaframe[i, 6:13] %in% areaframe$x_y[areaframe$area]) &
        areaframe$area[i]) {
      areaframe[i, "frame"] <- T
    }
  }
  
  for (i in seq_len(length(type))) {
    if (is.na(type[i])) {
    } else if (type[i] == "tiff") {
      tiff(
        filename = paste0(savepath, "/", savename, ".tiff"),
        width = width * 72, height = height * 72,
        family = family, compression = "zip"
      )
    } else if (type[i] == "pdf") {
      pdf(
        file = paste0(savepath, "/", savename, ".pdf"),
        width = width, height = height,
        family = family
      )
    } else if (type[i] == "jpg") {
      jpeg(
        filename = paste0(savepath, "/", savename, ".jpg"),
        width = width * 72, height = height * 72,
        family = family
      )
    } else if (type[i] == "png") {
      png(
        filename = paste0(savepath, "/", savename, ".png"),
        width = width * 72, height = height * 72,
        family = family
      )
    }
    
    if (lightmode) {
      lightmode()
      colorscale <- colorRampPalette(c("white", color))(1000)
    } else {
      darkmode()
      colorscale <- colorRampPalette(c("black", color))(1000)
    }
    showtext::showtext_auto()
    print(
      image(mse,
            formula = tic ~ x * y,
            smooth.image = smooth.image,
            superpose = superpose,
            normalize.image = normalize.image,
            contrast.enhance = contrast.enhance,
            xlab = xlab,
            ylab = ylab,
            colorscale = rgb(55, 55, 55, 200, maxColorValue = 255),
            xlim = c(xmin - xlength * 0.05, xmax + xlength * 0.05),
            ylim = c(ymin - ylength * 0.1, ymax + ylength * 0.1),
            asp = asp,
            subset = areaframe$frame
      )
    )
    print(
      image(mse,
            formula = tic ~ x * y,
            smooth.image = smooth.image,
            superpose = superpose,
            normalize.image = normalize.image,
            contrast.enhance = contrast.enhance,
            xlab = xlab,
            ylab = ylab,
            colorscale = colorscale,
            xlim = c(xmin - xlength * 0.05, xmax + xlength * 0.05),
            ylim = c(ymin - ylength * 0.1, ymax + ylength * 0.1),
            asp = asp,
            subset = !areaframe$frame,
            add = T
      )
    )
    
    if (is.na(type[i])) {
    } else {
      dev.off()
    }
    showtext::showtext_auto(FALSE)
  }
  
  gc(reset = TRUE)
  return("完成绘图")
}

#' 根据mz批量选区质谱成像图
#'
#' @param filename 文件路径
#' @param areafile 区域文件路径
#' @param mz mz值，为空时自动提取
#' @param savepath 保存路径
#' @param savename 保存名称
#' @param type 图片格式
#' @param width 图片宽度
#' @param height 图片长度
#' @param asp 图像长宽比
#' @param lightmode 成像模式
#' @param fun 处理函数
#' @param smooth.image 是否平滑处理
#' @param superpose 是否分面绘制
#' @param normalize.image 归一化模式
#' @param contrast.enhance 成像强对比
#' @param xlab x轴标签
#' @param ylab y轴标签
#' @param color 成像颜色
#' @param family 字体
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param expand 绘制区域线时候，定义凸包周围的额外空间的范围。默认为0.1，当所选区域很薄，比如只有2个像素点宽度的时候，可能线绘制不出来，就可以将expand值调大。
#' @param ... 见[Cardinal::image()]
#'
#' @export
imzmlareaimage2 <- function(filename,
                            areafile,
                            mz = NULL,
                            savepath = "./",
                            savename = "test",
                            type = "jpg",
                            width = 8,
                            height = 7,
                            asp = 1,
                            lightmode = F,
                            fun = mean,
                            smooth.image = "gaussian",
                            superpose = F,
                            normalize.image = "linear",
                            contrast.enhance = "none",
                            xlab = "",
                            ylab = "",
                            color = c(
                              "blue2", "cyan2", "yellow",
                              "brown1", "firebrick3"
                            ),
                            family = "sans",
                            mass.range = NULL,
                            resolution = 5,
                            units = "ppm",
                            ...) {
  print("质谱成像绘图开始")
  
  # 判断path是否存在，不存在进行创建
  if (!dir.exists(savepath)) {
    dir.create(path = savepath, recursive = T)
  }
  
  suppressMessages(library("Cardinal"))
  
  mse <- readMSIData(filename,
                     attach.only = TRUE,
                     mass.range = mass.range, resolution = resolution, units = units
  )
  samplearea <- readRDS(file = areafile)
  
  if (length(mz) == 0) {
    realmz <- mz(mse)
  } else {
    realmz <- mz
  }
  
  xmax <- max(coord(mse)$x)
  xmin <- min(coord(mse)$x)
  ymax <- max(coord(mse)$y)
  ymin <- min(coord(mse)$y)
  xlength <- xmax - xmin + 1
  ylength <- ymax - ymin + 1
  
  areaframe <- data.frame(
    x = coord(mse)$x,
    y = coord(mse)$y,
    x_y = paste0(coord(mse)$x, "_", coord(mse)$y),
    area = samplearea,
    frame = F,
    up = paste0(coord(mse)$x, "_", coord(mse)$y + 1),
    down = paste0(coord(mse)$x, "_", coord(mse)$y - 1),
    left = paste0(coord(mse)$x - 1, "_", coord(mse)$y),
    right = paste0(coord(mse)$x + 1, "_", coord(mse)$y),
    upleft = paste0(coord(mse)$x - 1, "_", coord(mse)$y + 1),
    upright = paste0(coord(mse)$x + 1, "_", coord(mse)$y + 1),
    downleft = paste0(coord(mse)$x - 1, "_", coord(mse)$y - 1),
    downright = paste0(coord(mse)$x + 1, "_", coord(mse)$y - 1)
  )
  
  
  for (i in 1:dim(areaframe)[1]) {
    if (!all(areaframe[i, 6:13] %in% areaframe$x_y[areaframe$area]) &
        areaframe$area[i]) {
      areaframe[i, "frame"] <- T
    }
  }
  areaframe <<- areaframe
  
  for (j in seq_len(length(realmz))) {
    print(j)
    for (i in seq_len(length(type))) {
      if (is.na(type[i])) {
      } else if (type[i] == "tiff") {
        tiff(
          filename = paste0(
            savepath, "/", savename, "-",
            format(realmz[j], nsmall = 5,trim = T), ".tiff"
          ),
          width = width * 72, height = height * 72,
          family = family, compression = "zip"
        )
      } else if (type[i] == "pdf") {
        pdf(
          file = paste0(
            savepath, "/", savename, "-",
            format(realmz[j], nsmall = 5,trim = T), ".pdf"
          ),
          width = width, height = height,
          family = family
        )
      } else if (type[i] == "jpg") {
        jpeg(
          filename = paste0(
            savepath, "/", savename, "-",
            format(realmz[j], nsmall = 5,trim = T), ".jpg"
          ),
          width = width * 72, height = height * 72,
          family = family
        )
      } else if (type[i] == "png") {
        png(
          filename = paste0(
            savepath, "/", savename, "-",
            format(realmz[j], nsmall = 5,trim = T), ".png"
          ),
          width = width * 72, height = height * 72,
          family = family
        )
      }
      
      if (lightmode) {
        lightmode()
        colorscale <- colorRampPalette(c("white", color))(1000)
      } else {
        darkmode()
        colorscale <- colorRampPalette(c("black", color))(1000)
      }
      showtext::showtext_auto()
      try({
        print(
          image(mse,
                mz = realmz[j],
                fun = mean,
                smooth.image = smooth.image,
                superpose = superpose,
                normalize.image = normalize.image,
                contrast.enhance = contrast.enhance,
                xlab = xlab,
                ylab = ylab,
                colorscale = rgb(55, 55, 55, 200, maxColorValue = 255),
                xlim = c(xmin - xlength * 0.05, xmax + xlength * 0.05),
                ylim = c(ymin - ylength * 0.1, ymax + ylength * 0.1),
                asp = asp,
                subset = areaframe$frame
          )
        )
        
        print(
          image(mse,
                mz = realmz[j],
                fun = mean,
                smooth.image = smooth.image,
                superpose = superpose,
                normalize.image = normalize.image,
                contrast.enhance = contrast.enhance,
                xlab = xlab,
                ylab = ylab,
                colorscale = colorscale,
                xlim = c(xmin - xlength * 0.05, xmax + xlength * 0.05),
                ylim = c(ymin - ylength * 0.1, ymax + ylength * 0.1),
                asp = asp,
                subset = !areaframe$frame,
                add = T
          )
        )
      })
      
      
      if (is.na(type[i])) {
      } else {
        dev.off()
      }
      showtext::showtext_auto(FALSE)
    }
  }
  
  return("完成绘图")
}


make_spot_vertices <- function (spot_positions, vertex_offsets) {
  spot_vertices <- merge(spot_positions, vertex_offsets)
  spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
  spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
  as.data.frame(spot_vertices)
}

select_spot_positions <- function (cdata, x = "col", y = "row", fill = "spatial.cluster", ...) {
  assertthat::assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
  if (is.character(fill) && length(fill) == 1) {
    spot_positions <- cdata[, c(x, y, fill)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")
  }
  else if (is.vector(fill) || is.factor(fill)) {
    assertthat::assert_that(nrow(cdata) == length(fill))
    spot_positions <- cdata[, c(x, y)]
    colnames(spot_positions) <- c("x.pos", "y.pos")
    spot_positions$fill <- fill
  }
  spot_positions$spot <- rownames(spot_positions)
  spot_positions
}

make_square_spots <- function (cdata, fill = "spatial.cluster", ...) {
  spot_positions <- select_spot_positions(cdata, fill = fill, ...)
  # spot_positions <<- spot_positions
  x.offset <- unique(spot_positions[order(spot_positions[,1]),1])
  x.offset <- x.offset[-1] - x.offset[-length(x.offset)]
  x.offset <- min(x.offset)
  y.offset <- unique(spot_positions[order(spot_positions[,2]),2])
  y.offset <- y.offset[-1] - y.offset[-length(y.offset)]
  y.offset <- min(y.offset)
  vertex_offsets <- data.frame(x.offset = c(0, x.offset, x.offset, 0), 
                               y.offset = c(0, 0, y.offset, y.offset))
  # vertex_offsets <- vertex_offsets * scale.factor
  make_spot_vertices(spot_positions, vertex_offsets)
}

normalize.image.linear <- function(x, range = c(0, 1)){
  
  if (all(is.na(x))){return(x)}
  oldmin <- min(x, na.rm = TRUE)
  oldmax <- max(x, na.rm = TRUE)
  min <- range[1]
  max <- range[2]
  ((x - oldmin) * (max - min)/(oldmax - oldmin)) + min
  
}

image.gaussian <- function (x, window = 3, ...) {
  if (all(is.na(x))){return(x)}
  r <- floor(window/2)
  sd <- window/4
  x.new <- .Call("C_gaussianFilter", x, r, sd, PACKAGE = "Cardinal")
  x.new <- max(x, na.rm = TRUE) * x.new/max(x.new, na.rm = TRUE)
  x.new
}

#' @export
imzmlareaimage3 <- function(filename,
                            areafile,
                            savepath = "./",
                            savename = "test",
                            type = c("png","pdf"),
                            width = 8,
                            height = 7,
                            asp = 1,
                            lightmode = F,
                            fun = mean,
                            smooth.image = "gaussian",
                            normalize.image = "linear",
                            contrast.enhance = "none",
                            xlab = "",
                            ylab = "",
                            colors = c("blue2", "cyan2", "yellow","brown1", "firebrick3"),
                            Area_colors = c("#EE7621","#B452CD","#e31a1c","#EEC900","#00EE00","#33a02c","#00688B","#8B1A1A","#CD5B45","#008B8B","#5D478B","#FF34B3","#00FF7F","#008B45","#7A67EE"),
                            family = "sans",
                            mass.range = NULL,
                            resolution = 5,
                            units = "ppm",
                            vague = F,
                            iteration = 2,
                            maxPixels = 1000000,
                            mapmz = F,
                            mz = NULL,
                            filtermz = NULL,
                            linetype = 1,
                            expand = 0.1,
                            ...) {
  print("质谱成像绘图开始")
  
  suppressMessages(library("ggplot2"))
  suppressMessages(library("magrittr"))
  suppressMessages(library("ggforce"))
  suppressMessages(library("concaveman"))
  suppressMessages(library("Cardinal"))
  
  # data
  mse <- readdata(filename=filename,
                  attach.only = TRUE,
                  mass.range = mass.range, 
                  resolution = resolution, 
                  units = units,
                  filtermz = filtermz)
  samplearea <- sapply(areafile, readRDS) %>% data.frame(check.names = F)
  colnames(samplearea) <- gsub(pattern = "-neg.rds|-pos.rds", replacement = "", x = basename(areafile))
  
  # zone label 
  subset_coor <- data.frame(Cardinal::coord(mse), samplearea,check.names = F) %>% 
    tidyr::pivot_longer(., cols = -(1:2), names_to = "Area") %>%
    dplyr::filter(., value)
  # print(head(subset_coor))
  
  subset_coor$y <- subset_coor$y*asp
  coord(mse)$y <- coord(mse)$y*asp
  
  # image.gaussian and normalize.image.linear
  # spe <- Cardinal::spectra(mse)
  if(vague){
    mse <- imzmlvague(mse = mse,
                      iteration = iteration,
                      maxPixels = maxPixels,
                      smooth.image = smooth.image)
    smooth.image  <- "none"
  } 
  
  if(smooth.image != "none"){
    spe <- spatialsmooth(mse,window = 3,smooth.image = smooth.image)
  }else{
    spe <- spectra(mse)
  }
  
  if (lightmode) {
    background <- "white"
  } else {
    background <- "black"
  }
  
  if(mapmz){
    if (length(mz) == 0) {
      realmz <- mz(mse)
    } else {
      realmz <- mz
    }
  }else{
    realmz <- "all"
  }
  
  suppressWarnings(library(foreach))
  suppressWarnings(library(doParallel))
  registerDoParallel(cores=5)
  foreach(mz_ = realmz) %dopar% {
  # for (mz_ in realmz){
    if(mz_ == "all"){
      fill <- "all"
      legendtitle <- ""
      tic <- pixelApply(mse, fun)
      data_x_y <- data.frame(Cardinal::coord(mse), tic)
    }else{
      mz_which <- which(mz_==mz(mse))
      fill <- format(as.numeric(mz_), nsmall = 5, trim=T)
      legendtitle <- paste0("-",fill)
      data_x_y <- data.frame(Cardinal::coord(mse), t(matrix(spe, nrow = nrow(mse))))[,c(1:2, 2+mz_which)]
    }
    dat <- data_x_y
    colnames(dat)[3] <-  fill
    
    # coordinate data
    vertices <- make_square_spots(dat, fill=fill, x="x", y="y", scale.factor = 1)
    
    # plot
    if (background=="black") col <- "white" else col <- "black"
    gradientn <- colorRampPalette(c(background, colors))(1000) 
    Area_color <- Area_colors[1:length(areafile)]
    names(Area_color) <- unique(subset_coor$Area)
    
    splot <- ggplot() + 
      geom_polygon(data = vertices, aes(x = x.vertex, y = y.vertex, group = spot, fill = fill), 
                   linewidth = 1/10^9) + scale_size(range = c(0, 1/10^8))+
      coord_equal() + theme_void(base_line_size = 1/10^9, base_rect_size = 1/10^9) +
      scale_fill_gradientn(colours = gradientn) + 
      scale_y_reverse() + 
      guides(fill=guide_colorbar(title = "intensity"))+
      theme(plot.margin=unit(rep(1,4),'cm'),
            legend.text = element_text(colour = col),
            legend.title =  element_text(colour = col)) +
      geom_mark_hull(data=subset_coor, aes(x=x+0.5, y=y+0.5, color=Area),  
                     linewidth=1, linetype=linetype, show.legend = T, con.cap = 0, 
                     expand = unit(expand, "mm"),radius = unit(1, "mm"),concavity = 3.5) +
      scale_color_manual(values = Area_color)+
      guides(color = guide_legend(order = 1))
    
    ggplotsave(plot = splot,
               savepath = savepath, 
               mapname = paste0(savename,legendtitle),
               imagetype = type,
               bg = background,
               family = family,
               width = width,
               height = height)
  }
  
  return("完成绘图")
}
