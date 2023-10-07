#!/opt/conda/bin/Rscript

#2.1 delete background
#' @export
delete_figbg <- function(file,
                         savepath=NULL,
                         ps = 5,
                         tv = 0.1,
                         bg_alpha = 0,
                         fig_alpha = 0.2,
                         ...
){
  if (!tools::file_ext(basename(file))=="png"){
    image <- magick::image_read(file)
    images_png <- magick::image_convert(image, "png")
    file <- paste0(tools::file_path_sans_ext(file), ".png")
    magick::image_write(images_png, path = file, format = "png")
  }

  x <- png::readPNG(source=file, ...)
  dimx <- dim(x)
  message(paste0("'",file,"' dimension:",
                 " height-",dimx[1],
                 ", width-", dimx[2],
                 ", channel-", dimx[3]
  ))

  if (dimx[3] %in% 3:4){

    n <- dimx[1]*dimx[2]
    r <- x[1:n]
    g <- x[(n+1):(2*n)]
    b <- x[(2*n+1):(3*n)]

    ps <- dimx[1]*(ps-1) + ps
    sel <- abs(r-r[ps])<tv & abs(g-g[ps])<tv & abs(b-b[ps])<tv
    alpha <- rep(1, n)
    alpha[sel] <- bg_alpha
    alpha[!sel] <- fig_alpha
    x[(3*n+1):(4*n)] <- alpha

    if (dimx[3]==3){
      x <- array(x, dim =c(dimx[1], dimx[2], 4) )
    }
  }else(stop("The channel's number must be 3 or 4, please check the image file."))

  file_name <- paste0(tools::file_path_sans_ext(basename(file)), "_debg.png")

  if(length(savepath)>0){
    if(!dir.exists(savepath)){dir.create(savepath,recursive=T)}
    file_name <- paste(savepath, file_name, sep = "/")
  }

  print("delete_figbg Done!")
  png::writePNG(x, file_name)
  return(file_name)
}

#2.2 cut off figure
#' @export
fig_cut <- function(file,
                    savepath=NULL,
                    pixel_thresh=0.85, num=5,
                    ymax_adj=0, ymin_adj=0,    # adjusting for coordinate
                    xmax_adj=0, xmin_adj=0,
                    interaction=T,
                    rotate_angle=0,
                    flip=F,
                    flop=F,
                    resize=F,
                    width=1000,
                    height=1000){
  # read image
  ebi <- EBImage::readImage(file)
  w <- dim(ebi)[1]
  h <- dim(ebi)[2]
  dims <- w*h

  # filter background
  ebi_thresh <- ebi>pixel_thresh
  data <- EBImage::imageData(ebi_thresh)
  pixel <- which(!data[1:dims])

  ## cut coordinate
  y_list <- c()
  x_list <- c()
  for(i in 1:length(pixel)){
    y <- floor(pixel[i]/w)
    x <- pixel[i]-w*y
    y_list[i] <- y
    x_list[i] <- x
  }

  y_num <- rle(sort(y_list))
  y_fine <- y_num$values[which(y_num$lengths>num)]
  ymax <- max(y_fine)+ymax_adj
  ymin <- min(y_fine)-ymin_adj

  x_num <- rle(sort(x_list))
  x_fine <- x_num$values[which(x_num$lengths>num)]
  xmax <- max(x_fine)+xmax_adj
  xmin <- min(x_fine)-xmin_adj

  # interaction for checking cutting location
  if (!is.null(dev.list())){
    dev.off()
  }

  if (interaction){
    imgjpg <- imager::load.image(file)
    plot(imgjpg,xlim = c(1,imager::width(imgjpg)),ylim = c(imager::height(imgjpg),1))
    rect(xleft = xmin, ybottom = ymax, xright = xmax, ytop = ymin, lwd=2,lty=1, border = "red")

    x<-readline(prompt = "Whether to shut down the graphic device (yes or no): ")
    if (x=="yes"){
      dev.off()
    }
  }else{
    imgjpg <- imager::load.image(file)
    filename <- paste(tools::file_path_sans_ext(basename(file)),pixel_thresh, num, "maker_cut.png", sep = "_")

    if(length(savepath)>0){
      if(!dir.exists(savepath)){dir.create(savepath,recursive=T)}
      filename <- paste(savepath, filename, sep = "/")
    }

    png(filename, width = imager::width(imgjpg), height = imager::height(imgjpg))
    plot(imgjpg,xlim = c(1,imager::width(imgjpg)),ylim = c(imager::height(imgjpg),1))
    rect(xleft = xmin, ybottom = ymax, xright = xmax, ytop = ymin, lwd=2,lty=1, border = "red")
    dev.off()
  }

  # cutting matrix and remodeling figure
  loc <- c()
  for (i in 1:dims){
    y <- floor(i/w)
    x <- i-w*y
    loc[i] <- ifelse((y>=ymin & y<=ymax)&(x>=xmin & x<=xmax),T,F)
  }

  data_array <- EBImage::imageData(ebi)
  frames.total <- dim(ebi)[3]
  r <- data_array[1:dims][loc]
  g <- data_array[(dims+1):(2*dims)][loc]
  b <- data_array[(2*dims+1):(3*dims)][loc]
  a <- data_array[(3*dims+1):(4*dims)][loc]
  if (frames.total==2){
    figure_cut <- array(c(r,g), dim=c(xmax-xmin+1,ymax-ymin+1, 2))
  }else if(frames.total==3){
    figure_cut <- array(c(r,g,b), dim=c(xmax-xmin+1,ymax-ymin+1, 3))
  }else{
    figure_cut <- array(c(r,g,b,a), dim=c(xmax-xmin+1,ymax-ymin+1, 4))
  }

  EBImage::imageData(ebi) <- figure_cut

  # rotate or flip and other operation
  if (flip) ebi <- EBImage::flip(ebi)
  if (flop) ebi <- EBImage::flop(ebi)
  if (rotate_angle!=0) ebi <- EBImage::rotate(ebi, rotate_angle, bg.col = "white")
  if (resize) ebi <- EBImage::resize(ebi, w=width, h=height)

  # save
  file_name <- paste(tools::file_path_sans_ext(basename(file)),pixel_thresh, num, "cut.png", sep = "_")
  if(length(savepath)>0){
    if(!dir.exists(savepath)){dir.create(savepath,recursive=T)}
    file_name <- paste(savepath, file_name, sep = "/")
  }

  EBImage::writeImage(ebi, file_name, type="png")

  # Binary, indicating if the spot falls inside (1) or outside (0) of tissue
  in_tissue <- sapply(!data[1:dims], FUN = function(x)return(ifelse(x, 1, 0)))[loc]

  result <- list(figcut=ebi,
                 in_tissue=in_tissue,
                 file_name=file_name,
                 cut_coordinate=data.frame(xleft = xmin,
                                           ybottom = ymax,
                                           xright = xmax,
                                           ytop = ymin)
  )

  message("fig_cut Done!")

  return(result)
}

#' @export
fig_deal <- function(file,
                     savepath = NULL,
                     ps = 5,
                     tv = 0.1,
                     bg_alpha = 0,
                     fig_alpha = 0.2,
                     pixel_thresh = 0.85, num = 5,
                     ymax_adj = 0, ymin_adj = 0,
                     xmax_adj = 0, xmin_adj = 0,
                     ymax_adj_ = 0, ymin_adj_ = 0,
                     xmax_adj_ = 0, xmin_adj_ = 0,
                     interaction = F,
                     rotate_angle = 0,
                     flip = F,
                     flop = F,
                     resize = F,
                     width = 1000,
                     height = 1000){
  library(dplyr)

  result <- delete_figbg(file,
                         savepath = savepath,
                         ps = ps,
                         tv = tv,
                         bg_alpha = bg_alpha,
                         fig_alpha = fig_alpha) %>% fig_cut(.,
                                                            savepath = savepath,
                                                            pixel_thresh=pixel_thresh, num=num,
                                                            ymax_adj=ymax_adj, ymin_adj=ymin_adj,
                                                            xmax_adj=xmax_adj, xmin_adj=xmin_adj,
                                                            interaction=interaction,
                                                            rotate_angle=rotate_angle,
                                                            flip=flip,
                                                            flop=flop,
                                                            resize=resize,
                                                            width=width,
                                                            height=height)

  if (rotate_angle!=0){
    result <- delete_figbg(result$file_name,
                           savepath = savepath,
                           ps = ps,
                           tv = tv,
                           bg_alpha = bg_alpha,
                           fig_alpha = fig_alpha) %>% fig_cut(.,
                                                              savepath = savepath,
                                                              interaction=interaction,
                                                              pixel_thresh=pixel_thresh, num=num,
                                                              ymax_adj=ymax_adj_, ymin_adj=ymin_adj_,
                                                              xmax_adj=xmax_adj_, xmin_adj=xmin_adj_)
  }

  return(result)
}

